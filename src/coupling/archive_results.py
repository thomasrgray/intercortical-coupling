"""Move previous generated results into a timestamped archive; preview by default.

Only known output folders and figure artifacts are selected. Raw recordings,
session metadata and source code are not selected. The default archive root
is outside the selected output trees.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

SESSION_OUTPUTS = ('LFP_channel_check', 'cleaned_LFP_artifacts', 'LFP_analyses')
FIGURE_SUFFIXES = {'.png', '.svg', '.pdf', '.jpg', '.jpeg', '.npz', '.npy', '.csv', '.json', '.pkl'}
PROTECTED_SUFFIXES = {'.dat', '.rhd', '.rhs', '.h5', '.hdf5', '.info', '.params'}


def session_paths(data_list=None, datasets=()):
    paths = list(datasets)
    if data_list:
        paths += [line.strip() for line in Path(data_list).expanduser().read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith('#')]
    result = sorted({Path(p).expanduser().resolve(strict=True) for p in paths})
    if any(not p.is_dir() for p in result):
        raise ValueError('Every session must be a directory')
    return result


def plan_archive(repo, sessions, extra_outputs=()):
    """Return (source, archive-relative destination) pairs, without changing files."""
    repo = Path(repo).resolve(strict=True)
    if not (repo/'pyproject.toml').is_file() or not (repo/'src/coupling').is_dir():
        raise ValueError('repo-root must be the intercortical-coupling checkout')
    pairs = []
    if (repo/'results').exists():
        pairs.append((repo/'results', Path('repository/results')))
    # Figure folders contain source READMEs as well as generated artifacts.
    # Archive only artifact files, preserving all code and documentation.
    figures = repo/'figures'
    if figures.exists():
        pairs.extend((p, Path('repository')/p.relative_to(repo)) for p in sorted(figures.rglob('*'))
                     if p.is_file() and p.suffix.lower() in FIGURE_SUFFIXES
                     and not any('archive' in part.lower() for part in p.relative_to(figures).parts))
    for session in sessions:
        identity = session.name + '_' + hashlib.sha256(str(session).encode()).hexdigest()[:10]
        pairs.extend((session/name, Path('sessions')/identity/name) for name in SESSION_OUTPUTS
                     if (session/name).exists())
    for item in extra_outputs:
        path = Path(item).expanduser().absolute()
        if not path.exists():
            raise ValueError(f'Output directory does not exist: {path}')
        # Explicit output options still cannot select a source tree or session root.
        if path.resolve() == repo or path.resolve() in sessions or (path/'.git').exists():
            raise ValueError(f'Not a generated output directory: {path}')
        identity = path.name + '_' + hashlib.sha256(str(path).encode()).hexdigest()[:10]
        pairs.append((path, Path('extra')/identity))
    # Collapse duplicates/children if a parent directory is already selected.
    selected = []
    for source, dest in sorted(pairs, key=lambda pair: (len(pair[0].parts), str(pair[0]))):
        if any(source == p or p in source.parents for p, _ in selected):
            continue
        selected.append((source, dest))
    for source, _ in selected:
        entries = [source] + (list(source.rglob('*')) if source.is_dir() else [])
        if any(p.is_symlink() for p in entries):
            raise ValueError(f'Refusing symlink in selected output: {source}')
        if any(p.suffix.lower() in PROTECTED_SUFFIXES or p.name in
               {'trial_info_frame.csv', 'dig_in_frame.csv', '.git'} for p in entries):
            raise ValueError(f'Recording/metadata found inside output; inspect before archiving: {source}')
        if any(p.suffix.lower() in {'.py', '.sh', '.toml'} for p in entries):
            raise ValueError(f'Source code found inside selected output: {source}')
    return selected


def archive(pairs, destination):
    destination = Path(destination).absolute()
    for source, _ in pairs:
        if source == destination or source in destination.parents or destination in source.parents:
            raise ValueError('Archive destination must be outside every selected output')
    destination.mkdir(parents=True, exist_ok=False)
    manifest = destination/'manifest.json'
    records = [{'source': str(source), 'archived': str(destination/relative), 'status': 'pending'}
               for source, relative in pairs]
    payload = {'created_utc': datetime.now(timezone.utc).isoformat(), 'moves': records}
    def save():
        temporary = manifest.with_suffix('.tmp')
        temporary.write_text(json.dumps(payload, indent=2))
        temporary.replace(manifest)
    save()
    for record in records:
        target = Path(record['archived'])
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            raise FileExistsError(target)
        record['status'] = 'moving'
        save()
        shutil.move(record['source'], target)
        record['status'] = 'moved'
        save()
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path.cwd())
    parser.add_argument('--data-list')
    parser.add_argument('--dataset', action='append', default=[])
    parser.add_argument('--output-dir', action='append', default=[], help='Additional generated output directory; repeatable')
    parser.add_argument('--archive-root', type=Path, help='Default: sibling intercortical-coupling_archives folder')
    parser.add_argument('--apply', action='store_true', help='Move files; without this option only print the plan')
    args = parser.parse_args(argv)
    repo = args.repo_root.expanduser().resolve()
    try:
        pairs = plan_archive(repo, session_paths(args.data_list, args.dataset), args.output_dir)
        root = (args.archive_root.expanduser().resolve() if args.archive_root else
                repo.parent/(repo.name+'_archives'))
        destination = root/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
        for source, relative in pairs:
            print(f'{source} -> {destination/relative}')
        if not pairs:
            print('No previous results found.')
        elif args.apply:
            print(f'Archive complete: {archive(pairs, destination)}')
        else:
            print(f'Preview: {len(pairs)} paths. Re-run with --apply to move them.')
    except (ValueError, OSError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
