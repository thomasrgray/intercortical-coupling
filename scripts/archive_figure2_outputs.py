"""Move previous Figure 2 result directories into a dated archive; never delete."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import shutil


def archive_outputs(run):
    run = Path(run).expanduser().resolve()
    if not run.is_dir():
        raise ValueError(f'Run directory not found: {run}')
    # Only names used by the Figure 2 commands. Never touch preprocessing caches
    # or new figure2_preferred_* runs. Move complete outputs, not just PNGs.
    candidates = sorted(p for p in run.iterdir() if p.is_dir() and not p.is_symlink()
                        and (p.name == 'fig02_recording_validation' or
                             p.name.startswith(('figure2_comparison_', 'figure2_overlays_', 'figure2_readable_'))))
    if not candidates:
        print('No previous Figure 2 output folders to archive.')
        return None
    destination = run/'archive'/('figure2_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    destination.mkdir(parents=True, exist_ok=False)
    manifest = dict(run=str(run), status='in_progress', folders=[])
    manifest_path = destination/'manifest.json'
    def save():
        manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    save()
    for source in candidates:
        files = [str(p.relative_to(source)) for p in sorted(source.rglob('*')) if p.is_file()]
        target = destination/source.name
        shutil.move(str(source), str(target))
        manifest['folders'].append(dict(original=str(source), archived=str(target), files=files))
        save()
        print(f'Archived {source.name} -> {target}')
    manifest['status'] = 'complete'
    save()
    return destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', help='Completed run containing prepared/ and respiration_audit/')
    archive_outputs(parser.parse_args().run)
