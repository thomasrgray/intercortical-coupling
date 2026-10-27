import json
from pathlib import Path
import pytest
from coupling.archive_results import archive, main, plan_archive


def fixture_paths(tmp_path):
    repo = tmp_path/'repo'
    (repo/'src/coupling').mkdir(parents=True)
    (repo/'pyproject.toml').touch()
    (repo/'results').mkdir()
    (repo/'results/old.npz').write_bytes(b'old analysis')
    (repo/'figures/panel').mkdir(parents=True)
    (repo/'figures/panel/README.md').write_text('keep')
    (repo/'figures/panel/plot.png').write_bytes(b'image')
    session = tmp_path/'TG22/session'
    (session/'LFP_analyses').mkdir(parents=True)
    (session/'LFP_analyses/old.svg').write_text('old plot')
    for name in ['recording.h5', 'board-ADC-00.dat', 'trial_info_frame.csv', 'session.info', 'session.params']:
        (session/name).write_bytes(b'original')
    return repo, session


def test_dry_run_then_archive_preserves_recordings_and_readmes(tmp_path, capsys):
    repo, session = fixture_paths(tmp_path)
    main(['--repo-root', str(repo), '--dataset', str(session)])
    assert 'Preview:' in capsys.readouterr().out
    assert (repo/'results/old.npz').exists()
    pairs = plan_archive(repo, [session])
    manifest = archive(pairs, tmp_path/'archive/run')
    records = json.loads(manifest.read_text())['moves']
    assert all(r['status'] == 'moved' and Path(r['archived']).exists() for r in records)
    assert not (repo/'results').exists()
    assert not (session/'LFP_analyses').exists()
    assert (repo/'figures/panel/README.md').read_text() == 'keep'
    assert (session/'recording.h5').read_bytes() == b'original'
    assert (session/'board-ADC-00.dat').read_bytes() == b'original'


def test_rejects_recordings_symlinks_and_nested_archive(tmp_path):
    repo, session = fixture_paths(tmp_path)
    pairs = plan_archive(repo, [session])
    with pytest.raises(ValueError, match='outside'):
        archive(pairs, repo/'results/archive')
    (repo/'results/raw.dat').touch()
    with pytest.raises(ValueError, match='Recording'):
        plan_archive(repo, [session])
    (repo/'results/raw.dat').unlink()
    (repo/'results/link').symlink_to(session, target_is_directory=True)
    with pytest.raises(ValueError, match='symlink'):
        plan_archive(repo, [session])


def test_collision_and_deduplication(tmp_path):
    repo, session = fixture_paths(tmp_path)
    pairs = plan_archive(repo, [session], [repo/'results'])
    assert len([s for s, _ in pairs if s == repo/'results']) == 1
    target = tmp_path/'existing'
    target.mkdir()
    with pytest.raises(FileExistsError):
        archive(pairs, target)
    assert (repo/'results/old.npz').exists()
