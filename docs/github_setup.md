# Create the remote

Prepared name: `intercortical-coupling`. Suggested owner: `thomasrgray`.
Keep the repository private while developing the manuscript. No remote was created
by this starter package; the available GitHub connector lacks that operation.

Create an empty private repository named `intercortical-coupling` on GitHub,
without initializing a README. Then, from the unpacked starter directory:

```bash
git init -b main
git add .
git commit -m "Initialize figure-first manuscript analysis"
git remote add origin https://github.com/thomasrgray/intercortical-coupling.git
git push -u origin main
```

Alternatively create it using an authenticated GitHub CLI after the first commit:
`gh repo create thomasrgray/intercortical-coupling --private --source=. --remote=origin --push`.
No raw recordings, derived arrays or synthetic results should be committed.
No open-source license has been selected on the lab's behalf.

## PyCharm route

1. Extract the ZIP and open the `intercortical-coupling` folder as a project.
2. Select a Python 3.10+ interpreter or create a virtual environment.
3. In PyCharm's terminal run `python -m pip install -e .`.
4. Run `python -m unittest discover -s tests -v` and the demo commands in README.
5. Enable Git through VCS / Enable Version Control Integration (wording varies).
6. Review the files in Commit: `data/`, `results/` and large recording files must
   remain excluded. Commit the source/docs/configs/tests.
7. Use the IDE's GitHub sharing option to create a PRIVATE repository, or add the
   empty repository URL as origin and push. No remote is included in this ZIP.
