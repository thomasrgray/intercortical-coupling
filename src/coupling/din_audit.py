"""Read-only pulse-width audit of every DIN across a curated session list."""
import argparse
import csv
import json
from pathlib import Path
import re
import numpy as np
from .events import read_binary, sample_rate

SUMMARY = ["session", "path", "filename", "din", "stimulus", "expected_route",
           "fs_hz", "recording_s", "complete_pulses", "partial_pulses",
           "min_ms", "p05_ms", "median_ms", "p95_ms", "max_ms",
           "width_class", "route_check", "status", "notes"]
PULSES = ["session", "filename", "din", "stimulus", "pulse_index",
          "on_sample", "off_sample_exclusive", "on_s", "observed_width_ms",
          "left_censored", "right_censored", "complete"]
ROUTES = {"oBA": "ortho", "oIA": "ortho", "BK": "ortho",
          "BA": "retro", "IA": "retro", "S": "retro"}


def read_pulses(path, chunk=1_000_000):
    """Report partial boundary pulses, never count their widths as complete."""
    if chunk < 1:
        raise ValueError("chunk must be positive")
    a = read_binary(path, "<u2")
    starts, stops, prev = [], [], 0
    for k in range(0, len(a), chunk):
        x = np.asarray(a[k:k+chunk])
        if np.any(x > 1):
            raise ValueError("Expected unpacked uint16 0/1 DIN; bitfields are unsupported")
        delta = np.diff(x.astype(np.int8), prepend=prev)
        starts.extend((np.flatnonzero(delta == 1)+k).tolist())
        stops.extend((np.flatnonzero(delta == -1)+k).tolist())
        prev = int(x[-1])
    if prev:
        stops.append(len(a))
    on, off = np.asarray(starts, dtype=np.int64), np.asarray(stops, dtype=np.int64)
    complete = (on > 0) & (off < len(a))
    return on, off, complete, len(a)


def labels_for(session):
    paths = sorted(session.glob("*.info"))
    if len(paths) != 1:
        return {}, f"Expected one .info; found {len(paths)}; labels unavailable", None
    try:
        info = json.loads(paths[0].read_text())
        taste = info["taste_params"]
        dins, labels = taste["dig_in_nums"], taste["tastes"]
        ids = [int(d) for d in dins]
        if len(ids) != len(labels) or len(ids) != len(set(ids)):
            raise ValueError("Ambiguous DIN/label lengths or duplicate DINs")
        if any(not isinstance(label, str) for label in labels):
            raise ValueError("Stimulus names must be strings")
        return dict(zip(ids, labels)), "", str(paths[0])
    except (ValueError, KeyError, TypeError, OSError) as exc:
        return {}, f"Invalid stimulus metadata: {exc}; labels unavailable", str(paths[0])


def session_paths(listing):
    listing = Path(listing).resolve()
    paths = [(listing.parent/line.strip()).resolve()
             for line in listing.read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    if not paths or len(paths) != len(set(paths)):
        raise ValueError("Session list must be nonempty and contain no duplicate paths")
    return paths


def classify(widths, label, split_ms):
    if not len(widths):
        return "", "no complete pulses"
    kind = "long" if np.median(widths) >= split_ms else "short"
    route = ROUTES.get(label)
    if route is None:
        return kind, "unassigned"
    expected = "long" if route == "ortho" else "short"
    return kind, "consistent" if kind == expected else "REVIEW"


def plot_results(rows, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    valid = [r for r in rows if r.get("median_ms", "") != ""]
    if not valid:
        return
    sessions = list(dict.fromkeys(r["path"] for r in rows))
    dins = sorted({r["din"] for r in valid})
    # Paginate so large session lists remain readable.
    for start in range(0, len(sessions), 15):
        page = sessions[start:start+15]
        values = np.full((len(page), len(dins)), np.nan)
        labels = {}
        for r in valid:
            if r["path"] in page:
                index = (page.index(r["path"]), dins.index(r["din"]))
                values[index] = r["median_ms"]
                labels[index] = f'{r["median_ms"]:.1f}\n{r["stimulus"] or "unassigned"}'
        fig, ax = plt.subplots(figsize=(max(7, len(dins)*1.1), max(3, len(page)*.65)))
        lo, hi = min(r["median_ms"] for r in valid), max(r["median_ms"] for r in valid)
        if lo == hi:
            lo, hi = lo/2, hi*2
        mesh = ax.imshow(np.ma.masked_invalid(values), aspect="auto",
                         norm=LogNorm(vmin=lo, vmax=hi), cmap="viridis")
        for (i, j), label in labels.items():
            ax.text(j, i, label, ha="center", va="center", fontsize=7,
                    color="black", bbox=dict(facecolor="white", alpha=.8, edgecolor="none"))
        ax.set_xticks(range(len(dins)), [f"DIN {d:02d}" for d in dins])
        ax.set_yticks(range(len(page)), [f"{start+i+1}: {Path(s).name}" for i, s in enumerate(page)])
        ax.set_title("Complete-pulse median duration; labels from .info")
        fig.colorbar(mesh, ax=ax, label="Pulse width (ms; logarithmic scale)")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(output/f"duration_comparison_{start//15+1:02d}.{ext}", dpi=180)
        plt.close(fig)


def summary_arrays(rows):
    """Keep every session, including sessions with no valid pulse measurements."""
    sessions = list(dict.fromkeys(r["path"] for r in rows))
    dins = sorted({int(r["din"]) for r in rows if r.get("din", "") != ""})
    widths = np.full((len(sessions), len(dins)), np.nan)
    counts = np.full_like(widths, np.nan)
    cells = {}
    for row in rows:
        if row.get("din", "") == "":
            continue
        key = (sessions.index(row["path"]), dins.index(int(row["din"])))
        # Duplicate physical-DIN error rows are not measurements.
        if row.get("median_ms", "") != "":
            widths[key] = float(row["median_ms"])
        if row.get("complete_pulses", "") != "":
            counts[key] = int(row["complete_pulses"])
        cells[key] = row
    return sessions, dins, widths, counts, cells


def plot_all_sessions(rows, output):
    """One unpaginated two-panel overview; missing values stay blank, not zero."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    sessions, dins, widths, counts, cells = summary_arrays(rows)
    if not sessions:
        return
    # A placeholder column lets an all-failed audit still display every session.
    if not dins:
        widths = np.full((len(sessions), 1), np.nan)
        counts = widths.copy()
    fig, axes = plt.subplots(1, 2, sharey=True,
        figsize=(max(12, 2.1*max(1, len(dins))+6), max(4, .7*len(sessions)+2)),
        layout="constrained")
    finite = widths[np.isfinite(widths) & (widths > 0)]
    lo, hi = (float(finite.min()), float(finite.max())) if finite.size else (1., 2000.)
    if lo == hi:
        lo, hi = lo/2, hi*2
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#eeeeee")
    images = [
        axes[0].imshow(np.ma.masked_invalid(widths), aspect="auto", cmap=cmap,
                       norm=LogNorm(vmin=lo, vmax=hi)),
        axes[1].imshow(np.ma.masked_invalid(counts), aspect="auto", cmap=cmap,
                       vmin=0, vmax=max(1., float(np.nanmax(counts))) if np.isfinite(counts).any() else 1.)
    ]
    for i in range(len(sessions)):
        for j in range(max(1, len(dins))):
            row = cells.get((i, j), {})
            label = row.get("stimulus") or "unassigned"
            partial = int(row.get("partial_pulses") or 0)
            flag = " *" if row.get("route_check") == "REVIEW" else ""
            status = row.get("status", "")
            if np.isfinite(widths[i, j]):
                left = f"{widths[i,j]:.1f} ms{flag}\n{label}"
            else:
                left = "ERROR" if status == "error" else "—"
            if np.isfinite(counts[i, j]):
                right = f"{int(counts[i,j])} complete\n{partial} partial"
            else:
                right = "ERROR" if status == "error" else "—"
            for ax, text in zip(axes, (left, right)):
                ax.text(j, i, text, ha="center", va="center", fontsize=7,
                        bbox=dict(facecolor="white", alpha=.85, edgecolor="none"))
    for ax, title in zip(axes, ("A  Median complete-pulse duration",
                                "B  Complete and boundary-partial pulse counts")):
        ax.set_title(title, fontsize=11)
        ax.set_xticks(range(max(1, len(dins))),
                      [f"DIN {d:02d}" for d in dins] or ["No readable DIN"])
        ax.set_xlabel("Digital input")
    names = []
    for i, session in enumerate(sessions):
        session_rows = [r for r in rows if r["path"] == session]
        errors = sum(r["status"] == "error" for r in session_rows)
        suffix = f" [errors: {errors}]" if errors else ""
        names.append(f"{i+1}: {Path(session).name}{suffix}")
    axes[0].set_yticks(range(len(sessions)), names, fontsize=8)
    axes[0].set_ylabel("Session")
    fig.colorbar(images[0], ax=axes[0], label="Median width (ms; log scale)", shrink=.8)
    fig.colorbar(images[1], ax=axes[1], label="Complete pulse count", shrink=.8)
    fig.suptitle("All-session DIN summary\n"
                 "* Duration/label mismatch; partial pulses excluded from medians; "
                 "gray = unavailable; labels from .info", fontsize=11)
    try:
        for ext in ("png", "pdf", "svg"):
            fig.savefig(Path(output)/f"all_sessions_summary.{ext}", dpi=180)
    finally:
        plt.close(fig)


def audit(listing, output, split_ms=100.):
    if not np.isfinite(split_ms) or split_ms <= 0:
        raise ValueError("split_ms must be positive and finite")
    sessions = session_paths(listing)
    output = Path(output).resolve()
    if any(output == s or s in output.parents for s in sessions):
        raise ValueError("Output must be outside source sessions")
    output.mkdir(parents=True, exist_ok=False)
    rows, sources = [], []
    with (output/"pulse_widths.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PULSES)
        writer.writeheader()
        for index, session in enumerate(sessions, 1):
            print(f"[{index}/{len(sessions)}] {session.name}", flush=True)
            labels, note, info = labels_for(session)
            sources.append(dict(session=str(session), info=info))
            try:
                fs = sample_rate(session/"info.rhd")
                files = sorted(session.glob("board-DIN-*.dat"))
                if not files:
                    raise ValueError("No board-DIN-*.dat files found")
            except (ValueError, OSError) as exc:
                rows.append(dict(session=session.name, path=str(session),
                                 status="error", notes=str(exc)))
                continue
            ids = [int(m.group(1)) for p in files
                   if (m := re.fullmatch(r"board-DIN-(\d+)\.dat", p.name))]
            duplicate = {d for d in ids if ids.count(d) > 1}
            session_rows = []
            lengths = set()
            for path in files:
                row = dict(session=session.name, path=str(session), filename=path.name,
                           fs_hz=fs, status="ok", notes=note)
                try:
                    match = re.fullmatch(r"board-DIN-(\d+)\.dat", path.name)
                    if not match:
                        raise ValueError("Noncanonical DIN filename; not interpreted")
                    din = int(match.group(1))
                    if din in duplicate:
                        raise ValueError(f"Duplicate physical DIN {din}")
                    label = labels.get(din, "")
                    row.update(din=din, stimulus=label, expected_route=ROUTES.get(label, ""))
                    on, off, complete, n = read_pulses(path)
                    lengths.add(n)
                    widths = (off[complete]-on[complete])*1000/fs
                    kind, check = classify(widths, label, split_ms)
                    row.update(recording_s=n/fs, complete_pulses=int(complete.sum()),
                               partial_pulses=int((~complete).sum()),
                               width_class=kind, route_check=check)
                    if len(widths):
                        row.update(dict(zip(["min_ms", "p05_ms", "median_ms", "p95_ms", "max_ms"],
                                            np.percentile(widths, [0, 5, 50, 95, 100]).tolist())))
                    if check == "REVIEW":
                        row["status"] = "warning"
                        row["notes"] += "; duration disagrees with configured long/short heuristic"
                    if (~complete).any():
                        row["status"] = "warning"
                        row["notes"] += "; boundary pulse(s): observed widths excluded from statistics"
                    for i, (a, b, good) in enumerate(zip(on, off, complete)):
                        writer.writerow(dict(session=str(session), filename=path.name, din=din,
                            stimulus=label, pulse_index=i+1, on_sample=int(a),
                            off_sample_exclusive=int(b), on_s=a/fs, observed_width_ms=(b-a)*1000/fs,
                            left_censored=bool(a == 0), right_censored=bool(b == n), complete=bool(good)))
                    print(f"  DIN {din:02d} {label or 'unassigned':10s} "
                          f"n={len(widths)} median={row.get('median_ms', 'NA')} ms {check}", flush=True)
                except (ValueError, OSError) as exc:
                    row.update(status="error", notes=str(exc))
                session_rows.append(row)
            for din in sorted(set(labels)-set(ids)):
                session_rows.append(dict(session=session.name, path=str(session), din=din,
                    stimulus=labels[din], status="error", notes="Declared DIN has no canonical file"))
            if len(lengths) > 1:
                for row in session_rows:
                    if row["status"] != "error":
                        row["status"] = "warning"
                    row["notes"] += "; DIN file sample lengths differ"
            if note:
                for row in session_rows:
                    if row["status"] == "ok":
                        row["status"] = "warning"
            rows.extend(session_rows)
    with (output/"din_duration_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY)
        writer.writeheader()
        writer.writerows(rows)
    (output/"audit.json").write_text(json.dumps(dict(
        session_list=str(Path(listing).resolve()), split_ms=split_ms, sources=sources,
        route_by_label=ROUTES, output_status="tables_complete",
        interpretation="Widths check long/short consistency, not odor identity or physical delivery duration.",
        partial_pulses="Boundary high segments are censored and excluded from width statistics."
    ), indent=2))
    plot_results(rows, output)
    plot_all_sessions(rows, output)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("session_list")
    parser.add_argument("output")
    parser.add_argument("--split-ms", type=float, default=100.,
                        help="Heuristic long/short boundary in ms (default 100); never assigns identities")
    args = parser.parse_args()
    rows = audit(args.session_list, args.output, args.split_ms)
    print(f"Saved DIN audit to {Path(args.output).resolve()}")
    if any(row["status"] == "error" for row in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
