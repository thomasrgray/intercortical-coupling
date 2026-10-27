"""Run the pinned Julia LFP workflow, recording stage status and live logs."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-list', default='configs/fin_inter_list.txt')
    parser.add_argument('--dataset', action='append', default=[])
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--stages', nargs='+', choices=['validate', 'clean', 'plot', 'average', 'compare', 'phase'],
                        default=['validate', 'clean', 'plot', 'average', 'compare'])
    parser.add_argument('--exclude-taste', action='append', default=[])
    parser.add_argument('--mad-threshold', type=float, default=3)
    parser.add_argument('--animal-map')
    parser.add_argument('--permutations', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--phase-config', help='Reviewed phase config with input_backend=julia')
    parser.add_argument('--phase-bootstrap', type=int, default=500)
    args = parser.parse_args(argv)
    if args.jobs < 1 or args.permutations < 1 or args.mad_threshold <= 0:
        parser.error('jobs, permutations and mad-threshold must be positive')
    if 'phase' in args.stages:
        if not args.phase_config:
            parser.error('--phase-config is required for the phase stage')
        phase_config = Path(args.phase_config).expanduser().resolve(strict=True)
        configuration = json.loads(phase_config.read_text())
        if configuration.get('defaults', {}).get('input_backend') != 'julia':
            parser.error('The phase stage requires defaults.input_backend=julia')
        if args.phase_bootstrap < 2:
            parser.error('--phase-bootstrap must be at least 2')
    paths = args.dataset or [line.strip() for line in Path(args.data_list).expanduser().read_text().splitlines()
                             if line.strip() and not line.lstrip().startswith('#')]
    sessions = [Path(p).expanduser().resolve(strict=True) for p in paths]
    if not sessions or any(not p.is_dir() for p in sessions):
        parser.error('Provide at least one existing session directory')
    if len({p.name for p in sessions}) != len(sessions):
        parser.error('Session basenames must be unique')
    output = Path(args.output_dir).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('Use a new/empty output directory; archive previous results first')
    output.mkdir(parents=True, exist_ok=True)
    session_list = output/'sessions.txt'
    session_list.write_text(''.join(str(p)+'\n' for p in sessions))
    data_args = [arg for p in sessions for arg in ('--dataset', str(p))]
    steps = {
        'validate': ('validate_LFP_data', []),
        'clean': ('remove_lfp_artifacts', ['--mad-threshold', str(args.mad_threshold)] +
                  [arg for item in args.exclude_taste for arg in ('--exclude-taste', item)]),
        'plot': ('plot_LFP_spectrogram', ['--output-dir', str(output/'session_reports')]),
        'average': ('grand_average_LFP', ['--output-dir', str(output/'grand_average')]),
        'compare': ('compare_LFP_spectrograms', ['--output-dir', str(output/'statistics'),
                    '--permutations', str(args.permutations), '--seed', str(args.seed)] +
                    (['--animal-map', str(Path(args.animal_map).expanduser().resolve())] if args.animal_map else [])),
    }
    steps['phase'] = ('phase_coherence', [])
    report = {'source_commit': '1a4a4eb6346b00226e907228f4a2c3f658223907',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'sessions': list(map(str, sessions)), 'settings': vars(args), 'stages': []}
    manifest = output/'run.json'
    env = dict(os.environ, COUPLING_LFP_JOBS=str(args.jobs), MPLBACKEND='Agg', PYTHONUNBUFFERED='1')
    # Execute in canonical dependency order, regardless of option ordering.
    selected = [stage for stage in steps if stage in args.stages]
    for index, stage in enumerate(selected, 1):
        module, extra = steps[stage]
        command = [sys.executable, '-m', f'coupling.julia_lfp.{module}', *data_args, *extra]
        if stage == 'phase':
            command = [sys.executable, '-m', 'coupling.phase_coherence', str(session_list),
                       str(phase_config), str(output/'phase'), '--bootstrap', str(args.phase_bootstrap),
                       '--seed', str(args.seed)]
        status = {'stage': stage, 'command': command, 'status': 'running'}
        report['stages'].append(status)
        manifest.write_text(json.dumps(report, indent=2))
        print(f'[{index}/{len(selected)}] {stage}', flush=True)
        with (output/f'{stage}.log').open('w') as log:
            with subprocess.Popen(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, bufsize=1) as proc:
                try:
                    for line in proc.stdout:
                        print(line, end='', flush=True)
                        log.write(line)
                        log.flush()
                    code = proc.wait()
                except BaseException:
                    proc.terminate()
                    proc.wait()
                    status['status'] = 'interrupted'
                    manifest.write_text(json.dumps(report, indent=2))
                    raise
        status.update(status='complete' if code == 0 else 'failed', returncode=code)
        manifest.write_text(json.dumps(report, indent=2))
        if code:
            raise SystemExit(f'{stage} exited {code}; inspect {output / (stage + ".log")}. '
                             'Later stages were not run. Statistical stages also report insufficient animals as nonzero.')
    print(f'Completed. Run record: {manifest}')


if __name__ == '__main__':
    main()
