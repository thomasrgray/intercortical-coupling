"""Discover, validate and process multiple explicitly mapped Intan sessions."""
import argparse
import contextlib
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import traceback

from .events import inventory, sample_rate, select_six, write_csv
from .prepare import check_time, prepare
from .power import STIMULI


def write_json(path, value):
    path=Path(path)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2))
    tmp.replace(path)


def discover(root, destination):
    """Create an editable manifest; never infer animal identity from directory order."""
    root=Path(root).resolve();destination=Path(destination)
    if destination.exists(): raise FileExistsError(destination)
    if not root.is_dir(): raise ValueError(f'Not a session root: {root}')
    folders=sorted({p.parent for p in root.rglob('*.info')})
    sessions=[]
    for folder in folders:
        candidates=set(re.findall(r'(?:^|_)(TG\d+)(?=_|$)',folder.name))
        if re.fullmatch(r'TG\d+',folder.parent.name): candidates.add(folder.parent.name)
        animal=next(iter(candidates)) if len(candidates)==1 else ''
        sessions.append(dict(path=str(folder),session_id=folder.name,animal_id=animal))
    if not sessions: raise ValueError('No session .info files found (backups are ignored)')
    write_json(destination,dict(version=1,defaults=dict(reference_method='mean_car',notch_60hz=True,
        max_channels=5,min_trials=5,orders=[2,5,10,15,20],matched_pairs=[['BA','oBA'],['IA','oIA']]),sessions=sessions))
    return sessions


def resolve_file(folder, spec, key, pattern, errors, required=True):
    if key in spec:
        path=(folder/spec[key]).resolve()
        if path.parent != folder: errors.append(f'{key}: must be directly in session folder');return None
        matches=[path] if path.is_file() else []
    else: matches=sorted(folder.glob(pattern))
    if len(matches)!=1:
        if required or matches: errors.append(f'{key}: expected one {pattern}; found {len(matches)}. Set {key} explicitly if ambiguous.')
        return None
    return matches[0]


def flatten(values):
    for value in values:
        if isinstance(value,list): yield from flatten(value)
        else: yield int(value)


def preflight(spec, defaults):
    """Read metadata and scan every DIN/timestamp; report missing channels together."""
    folder=Path(spec['path']).resolve()
    report=dict(session_id=spec.get('session_id',''),animal_id=spec.get('animal_id',''),
        path=str(folder),errors=[],warnings=[],facts={},config=None)
    errors,warnings,facts=report['errors'],report['warnings'],report['facts']
    options={**defaults,**spec.get('options',{})}
    unknown=set(options)-{'reference_method','notch_60hz','max_channels','min_trials','orders','matched_pairs'}
    if unknown: errors.append('Unknown analysis options: '+', '.join(sorted(unknown)))
    if not report['animal_id']: errors.append('animal_id is required; set it in the manifest')
    if not report['session_id']: errors.append('session_id is required')
    if not folder.is_dir(): errors.append('Session folder does not exist');return report
    paths={key:resolve_file(folder,spec,key,pattern,errors,required) for key,pattern,required in [
        ('info','*.info',True),('layout','*electrode_layout*.csv',True),('params','*.params',False),
        ('header','info.rhd',True),('timestamps','time.dat',True)]}
    loaded={}
    for key in ('info','params'):
        if paths[key]:
            try:
                loaded[key]=json.loads(paths[key].read_text())
                if not isinstance(loaded[key],dict): raise ValueError('Expected JSON object')
            except (ValueError,OSError) as exc: errors.append(f'{key}: {exc}')
    if not paths['params']: warnings.append('No .params file: optional sampling-rate cross-check unavailable')
    info=loaded.get('info',{});params=loaded.get('params',{})
    if info.get('file_type') not in (None,'one file per channel'):
        errors.append('Only one file per channel Intan recordings are supported')
    fs=None
    if paths['header']:
        try:
            fs=sample_rate(paths['header']);facts['raw_fs_hz']=fs
            if fs!=int(fs): errors.append('Integer raw sampling rate required')
        except (ValueError,OSError) as exc: errors.append(f'header: {exc}')
    for label,value in [('info',info.get('recording_params',{}).get('sampling_rate')),('params',params.get('sampling_rate'))]:
        if value is not None and fs is not None and value!=fs: errors.append(f'{label} sampling_rate {value} disagrees with RHD {fs}')
    if params:
        facts['ignored_spike_bandpass_hz']=[params.get('bandpass_lower_cutoff'),params.get('bandpass_upper_cutoff')]
        warnings.append('Spike filters, clustering and auto_CAR flags in .params are not applied to raw LFP DAT')
    facts['acquisition_filters']=info.get('recording_params',{})
    # Anatomical mapping is explicit in .info; compare it with original CSV groups.
    region_map={};channels=[]
    try:
        for region,groups in info.get('electrode_layout',{}).items():
            for electrode in flatten(groups):
                if electrode in region_map: raise ValueError(f'Duplicate electrode {electrode} in .info')
                if region.upper() not in ('GC','PC'): raise ValueError(f'Unsupported region {region}')
                region_map[electrode]=region.upper()
        if not region_map: raise ValueError('Missing electrode_layout in .info')
        if paths['layout']:
            with paths['layout'].open() as h: table=list(csv.DictReader(h))
            for row in table:
                eid=int(row['electrode_ind']);region=region_map.get(eid)
                if region is None: raise ValueError(f'CSV electrode {eid} missing from .info')
                declared=(row.get('original_CAR_group') or row.get('CAR_group','')).upper()
                if declared in ('GC','PC') and declared!=region: raise ValueError(f'CSV/.info region disagreement for electrode {eid}')
                if row.get('original_CAR_group') and declared not in ('GC','PC'):
                    raise ValueError(f'Unrecognized original_CAR_group for electrode {eid}: {declared}')
                group=row.get('CAR_group','').strip()
                if not group: raise ValueError(f'Missing CAR_group for electrode {eid}')
                channels.append(dict(filename=row['filename'],electrode_id=eid,region=region,reference_group=group))
            if len({c['electrode_id'] for c in channels})!=len(channels): raise ValueError('Duplicate CSV electrode IDs')
            if len({c['filename'] for c in channels})!=len(channels): raise ValueError('Duplicate CSV amplifier filenames')
            if {c['electrode_id'] for c in channels}!=set(region_map): raise ValueError('CSV/.info electrode sets differ')
            if {c['region'] for c in channels}!={'GC','PC'}: raise ValueError('Need both GC and PC')
            groups={}
            for c in channels: groups.setdefault(c['reference_group'],[]).append(c)
            for group,members in groups.items():
                if len({c['region'] for c in members})!=1: raise ValueError(f'Reference group {group} spans regions')
                if options.get('reference_method','mean_car')=='mean_car' and len(members)<3:
                    raise ValueError(f'Reference group {group} has fewer than 3 channels before QC')
    except (KeyError,TypeError,ValueError,OSError) as exc: errors.append(f'electrode mapping: {exc}')
    facts['channels_by_region']={r:sum(c['region']==r for c in channels) for r in ('GC','PC')}
    for region,count in facts['channels_by_region'].items():
        if count<4: errors.append(f'{region}: need at least 4 mapped channels for the existing 3-peer QC rule')
    missing=[]
    for c in channels:
        path=(folder/c['filename']).resolve()
        if path.parent!=folder or not re.fullmatch(r'amp-[A-Za-z]-\d+\.dat',path.name):
            errors.append(f'Invalid amplifier filename: {c["filename"]}');continue
        if not path.is_file(): missing.append(c['filename'])
    if missing: errors.append('Missing mapped amplifier files: '+', '.join(missing))
    extras=sorted(p.name for p in folder.glob('amp-*.dat') if p.name not in {c['filename'] for c in channels})
    if extras: warnings.append('Unmapped amplifier files ignored: '+', '.join(extras))
    labels={};expected={}
    try:
        taste=info['taste_params'];dins=taste['dig_in_nums'];names=taste['tastes'];counts=taste['trial_count']
        if len(dins)!=6 or len(set(dins))!=6 or len(names)!=6 or set(names)!=set(STIMULI) or len(counts)!=6:
            raise ValueError('Need six unique DINs, the six supported stimulus names and six trial counts')
        labels={str(int(d)):s for d,s in zip(dins,names)}
        expected={int(d):int(n) for d,n in zip(dins,counts)}
        facts['stimulus_by_din']=labels;facts['expected_trials_by_din']=expected
    except (KeyError,ValueError,TypeError) as exc: errors.append(f'stimulus metadata: {exc}')
    if info.get('laser_params',{}).get('dig_in_nums'):
        errors.append('Optogenetic trials declared: batch pooling is blocked until an explicit opto-stratification policy is implemented')
    if fs is not None:
        try:
            records,rows=inventory(folder,fs);facts['din_inventory']=rows
            selected=select_six({d:len(v[0]) for d,v in records.items()})
            if set(selected)!=set(expected): errors.append('Measured six-DIN group disagrees with .info stimulus DINs')
            for din,count in expected.items():
                if din not in records or len(records[din][0])!=count: errors.append(f'DIN {din}: actual trial count disagrees with metadata ({count})')
            n=next(iter(records.values()))[2];facts['sample_count']=n;facts['duration_s']=n/fs
            facts['boundary_events']=sum(int(on<2*fs or on+5*fs>n) for din in selected for on in records[din][0])
            if facts['boundary_events']: warnings.append('Boundary events will be retained by ID but rejected by epoch QC')
            if paths['timestamps']: facts['first_timestamp']=check_time(paths['timestamps'],n)
            for c in channels:
                path=folder/c['filename']
                if path.is_file() and path.stat().st_size!=2*n: errors.append(f'Amplifier/DIN length mismatch: {path.name}')
            metadata_dins=info.get('dig_ins',{})
            for din,count in zip(metadata_dins.get('nums',[]),metadata_dins.get('trial_counts',[])):
                if din not in records or len(records[din][0])!=count: errors.append(f'Inventory DIN {din}: disagrees with .info dig_ins')
            other=[d for d in records if d not in selected and len(records[d][0])]
            if other: warnings.append(f'Non-stimulus DINs {other} are inventoried, not interpreted as opto or delivery events')
            facts['estimated_working_disk_gb']=len(channels)*(n/fs*1000*4+sum(expected.values())*7000*8)/1e9
        except (ValueError,OSError) as exc: errors.append(f'binary inventory/timestamps: {exc}')
    else:
        warnings.append('Binary timing checks require a valid info.rhd header')
    bad_dins=[p.name for p in folder.glob('board-DIN-*.dat') if not re.fullmatch(r'board-DIN-\d+\.dat',p.name)]
    if bad_dins: errors.append('Noncanonical DIN filenames (use original recording names): '+', '.join(bad_dins))
    if options.get('reference_method','mean_car') not in ('mean_car','none'): errors.append('reference_method must be mean_car or none')
    if not isinstance(options.get('notch_60hz',True),bool): errors.append('notch_60hz must be boolean')
    for key,minimum in [('max_channels',1),('min_trials',3)]:
        if not isinstance(options.get(key,minimum),int) or options.get(key,minimum)<minimum: errors.append(f'{key} must be an integer >= {minimum}')
    orders=options.get('orders',[2,5,10,15,20])
    if not isinstance(orders,list) or not orders or any(not isinstance(p,int) or p<1 or p>=125 for p in orders): errors.append('orders must be positive integers below 125 for 1.5 s VAR windows')
    pairs=options.get('matched_pairs',[['BA','oBA'],['IA','oIA']])
    if not isinstance(pairs,list) or not pairs or any(not isinstance(p,list) or len(p)!=2 or p[0] not in STIMULI or p[1] not in STIMULI or p[0]==p[1] for p in pairs): errors.append('Invalid matched_pairs; use [retro, ortho] stimulus names')
    facts['metadata_files']={k:str(p) if p else None for k,p in paths.items()}
    facts['metadata_sha256']={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items() if p and k in ('info','params','layout')}
    if not errors:
        report['config']=dict(animal_id=report['animal_id'],session_id=report['session_id'],raw_fs_hz=fs,
            header=paths['header'].name,timestamps=paths['timestamps'].name,channels=channels,
            stimulus_by_din=labels,reference_method=options.get('reference_method','mean_car'),
            notch_60hz=options.get('notch_60hz',True),batch_metadata=facts,**{k:options[k] for k in ('orders','matched_pairs') if k in options})
    report['options']=options
    return report


def load_manifest(path):
    path=Path(path).resolve();manifest=json.loads(path.read_text())
    if manifest.get('version')!=1: raise ValueError('Require manifest version 1')
    sessions=manifest.get('sessions',[])
    if not sessions: raise ValueError('Manifest has no sessions')
    ids=[];paths=[]
    for spec in sessions:
        spec['path']=str((path.parent/spec['path']).resolve())
        sid=spec.get('session_id','')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',sid): raise ValueError('session_id must be a safe, nonempty folder name')
        ids.append(sid);paths.append(spec['path'])
    if len(set(ids))!=len(ids) or len(set(paths))!=len(paths): raise ValueError('Duplicate session IDs or source folders')
    return manifest


def run_batch(manifest_path, output, execute=False):
    manifest=load_manifest(manifest_path);out=Path(output).resolve()
    for spec in manifest['sessions']:
        source=Path(spec['path'])
        if out==source or source in out.parents: raise ValueError('Batch output must be outside every raw session')
    if out.exists(): raise FileExistsError(out)
    out.mkdir(parents=True);write_json(out/'manifest_resolved.json',manifest)
    reports=[];summary=[]
    # Inspect all sessions before starting the first expensive preparation.
    for spec in manifest['sessions']:
        print('Checking '+spec['session_id'],flush=True)
        try: report=preflight(spec,manifest.get('defaults',{}))
        except Exception as exc:
            report=dict(session_id=spec['session_id'],animal_id=spec.get('animal_id',''),path=spec['path'],
                errors=[f'Malformed metadata or read error: {exc}'],warnings=[],facts={},config=None)
        reports.append(report)
        write_json(out/'preflight.json',reports)
    for report in reports:
        sid=report['session_id'];row=dict(session_id=sid,animal_id=report['animal_id'],status='blocked' if report['errors'] else 'ready',error='; '.join(report['errors']),output='',accepted_channel_events='',granger_fits_passed='',granger_fits_attempted='')
        work=out/sid;work.mkdir()
        if report['config']: write_json(work/'config.json',report['config'])
        if execute and not report['errors']:
            try:
                if shutil.disk_usage(out).free < report['facts'].get('estimated_working_disk_gb',0)*1e9*1.2:
                    raise ValueError('Insufficient disk for estimated working arrays plus 20% margin')
                row['status']='running';write_json(work/'status.json',row)
                with (work/'run.log').open('w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                    clean=prepare(report['path'],work/'config.json',work/'prepared')
                    from .lfp_coupling import analyze
                    options=report['options']
                    analyze(clean,work/'fig06_lfp_coupling',max_channels=options.get('max_channels',5),
                        min_trials=options.get('min_trials',5),matched_pairs=options.get('matched_pairs'),
                        orders=tuple(options.get('orders',[2,5,10,15,20])))
                with (work/'fig06_lfp_coupling'/'retention.csv').open() as h:
                    row['accepted_channel_events']=sum(int(r['accepted_trials']) for r in csv.DictReader(h))
                with (work/'fig06_lfp_coupling'/'granger_diagnostics.csv').open() as h:
                    diagnostics=list(csv.DictReader(h))
                row['granger_fits_passed']=sum(r['passed']=='True' for r in diagnostics)
                row['granger_fits_attempted']=len(diagnostics)
                row.update(status='completed',output=str(work/'fig06_lfp_coupling'))
            except Exception as exc:
                row.update(status='failed',error=str(exc))
                (work/'failure.txt').write_text(traceback.format_exc())
        summary.append(row);write_json(work/'status.json',row);write_csv(out/'summary.csv',summary)
        print(sid+': '+row['status'],flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    init=sub.add_parser('init',help='Discover session folders and write editable manifest')
    init.add_argument('root');init.add_argument('manifest')
    for name in ('check','run'):
        command=sub.add_parser(name);command.add_argument('manifest');command.add_argument('output')
    a=p.parse_args()
    if a.command=='init':
        sessions=discover(a.root,a.manifest);print(f'Wrote {len(sessions)} sessions; verify animal IDs and options in {a.manifest}')
    else:
        summary=run_batch(a.manifest,a.output,execute=a.command=='run')
        raise SystemExit(1 if any(r['status'] in ('failed','blocked') for r in summary) else 0)

if __name__=='__main__': main()
