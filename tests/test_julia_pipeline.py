import json
import sys
from pathlib import Path
from coupling.julia_lfp.pipeline import main


def test_phase_stage_uses_session_list_and_reviewed_config(tmp_path, monkeypatch):
    session = tmp_path/'TG22_session'
    session.mkdir()
    config = tmp_path/'phase.json'
    config.write_text(json.dumps({'defaults': {'input_backend': 'julia'}}))
    commands = []
    class Proc:
        stdout = iter(['done\n'])
        def __init__(self, cmd, **kwargs): commands.append(cmd)
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def wait(self): return 0
    monkeypatch.setattr('coupling.julia_lfp.pipeline.subprocess.Popen', Proc)
    out = tmp_path/'run'
    main(['--dataset',str(session),'--output-dir',str(out),'--stages','phase',
          '--phase-config',str(config),'--phase-bootstrap','20'])
    assert commands[0][:3] == [sys.executable, '-m', 'coupling.phase_coherence']
    assert Path(commands[0][3]).read_text().strip() == str(session)
    assert json.loads((out/'run.json').read_text())['stages'][0]['status'] == 'complete'
