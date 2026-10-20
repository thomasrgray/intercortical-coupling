"""Summarize saved LFP power and phase coupling across sessions and animals."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .figure6 import STYLE, COLORS, grid, save

COMPARISONS = {"route": ("ortho", "retro"), "BA": ("oBA", "BA"), "IA": ("oIA", "IA")}


def mean_finite(x, axis=0):
    x = np.asarray(x)
    n = np.isfinite(x).sum(axis=axis)
    return np.divide(np.nansum(x, axis=axis), n,
                     out=np.full(np.shape(n), np.nan), where=n > 0)


def animal_average(values, animals):
    """Sessions equally within animal, then animals equally; no pooled trials."""
    names = sorted(set(animals))
    if not names or any(not a or a.lower() in ("none", "unknown") for a in names):
        raise ValueError("Explicit animal IDs are required")
    per_animal = np.stack([mean_finite(values[np.array(animals) == a]) for a in names])
    return names, per_animal, mean_finite(per_animal)


def same_grid(a, b, keys):
    for key in keys:
        if not np.array_equal(a[key], b[key]):
            raise ValueError(f"Incompatible {key}; summarize matching runs separately")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_phase(run, comparison, sources, excluded):
    statuses = json.loads((run/"status.json").read_text())
    sources.append(run/"status.json")
    records, seen = [], set()
    for row in statuses:
        if row["status"] != "completed":
            excluded.append(dict(kind="phase", session=row["session"],
                                 reason=row.get("error") or row["status"]))
            continue
        if row["session"] in seen:
            raise ValueError("Duplicate session in phase status")
        seen.add(row["session"])
        folder = (run/row["output_key"]).resolve()
        if run not in folder.parents:
            raise ValueError("Phase output_key escapes run")
        provenance_path = folder/"provenance.json"
        provenance = json.loads(provenance_path.read_text())
        conditions_path = folder/"conditions.json"
        conditions = json.loads(conditions_path.read_text())
        pair = COMPARISONS[comparison]
        if any(conditions.get(label, {}).get("status") != "completed" for label in pair):
            excluded.append(dict(kind="phase", session=row["session"], reason="Missing completed matched conditions"))
            continue
        sources.extend([provenance_path, conditions_path])
        means, first = [], None
        for label in pair:
            path = (folder/conditions[label]["directory"]/"phase_coherence.npz").resolve()
            if folder not in path.parents:
                raise ValueError("Condition path escapes session output")
            with np.load(path, allow_pickle=False) as d:
                meta = {k: d[k] for k in ("time_s", "frequency_hz", "band", "band_limits_hz",
                                          "alignment_label", "selected_electrode_id")}
                grid(meta)
                draws = d["bootstrap_band"]
                if draws.ndim != 3 or draws.shape[0] < 2 or draws.shape[1:] != (len(meta["band"]), len(meta["time_s"])):
                    raise ValueError(f"Invalid bootstrap dimensions: {path}")
                if not np.isfinite(draws).all() or np.any(draws < 0) or np.any(draws > 1+1e-12):
                    raise ValueError(f"Invalid phase coherence: {path}")
                means.append(draws.mean(axis=0))
                del draws
            if first is not None:
                same_grid(first, meta, tuple(meta))
            else:
                first = meta
            sources.append(path)
        animal = provenance.get("animal_id") or row.get("animal_id")
        if not animal:
            raise ValueError(f"Missing animal ID for {row['session']}")
        signature = json.dumps({"profile": provenance["profile"], "settings": provenance["settings"]}, sort_keys=True)
        record = dict(session=Path(row["session"]).name, source_session=row["session"], animal=str(animal),
                      values=np.stack(means), meta=first, signature=signature)
        if records:
            same_grid(records[0]["meta"], first,
                      ("time_s", "frequency_hz", "band", "band_limits_hz", "alignment_label"))
            if signature != records[0]["signature"]:
                raise ValueError("Phase settings/profiles differ; summarize separately")
        records.append(record)
    return records


def read_power(root, comparison, sources):
    records, seen = [], set()
    for path in sorted(root.rglob("session_power.npz")):
        provenance_path = path.parent/"provenance.json"
        provenance = json.loads(provenance_path.read_text())
        with np.load(path, allow_pickle=False) as d:
            meta = {k: d[k] for k in ("time_s", "frequency_hz")}
            grid(meta)
            regions, stimuli = d["region"].tolist(), d["stimulus"].tolist()
            if len(set(regions)) != len(regions) or set(regions) != {"GC", "PC"} or len(set(stimuli)) != len(stimuli):
                raise ValueError(f"Invalid region/stimulus labels: {path}")
            x = d["z_power"]
            if x.shape != (len(regions), len(stimuli), len(meta["frequency_hz"]), len(meta["time_s"])) or np.isinf(x).any():
                raise ValueError(f"Invalid power axes/values: {path}")
            members = (("oBA", "oIA"), ("BA", "IA")) if comparison == "route" else tuple((s,) for s in COMPARISONS[comparison])
            # Ordinary mean requires both matched odors at each bin; never substitute one odor.
            values = np.stack([np.stack([
                np.mean(x[regions.index(region), [stimuli.index(s) for s in labels]], axis=0)
                for region in ("GC", "PC")]) for labels in members])
            if any(not np.isfinite(panel).any() for route in values for panel in route):
                raise ValueError(f"No finite power for a required route/region: {path}")
            animal, session = str(d["animal_id"].item()), str(d["session_id"].item())
            settings = json.loads(str(d["settings_json"].item()))
        if (animal, session) in seen:
            raise ValueError(f"Duplicate power session {animal}/{session}; use one run root")
        seen.add((animal, session))
        meta["alignment_label"] = f'Time from {provenance["preprocessing"]["event_alignment"]} (s)'
        signature = json.dumps(settings, sort_keys=True)
        if records:
            same_grid(records[0]["meta"], meta, ("time_s", "frequency_hz", "alignment_label"))
            if signature != records[0]["signature"]:
                raise ValueError("Power normalization/settings differ; summarize separately")
        records.append(dict(session=session, animal=animal, values=values, meta=meta, signature=signature))
        sources.extend([path, provenance_path])
    return records


def plot_power(records, out, comparison):
    values = np.stack([r["values"] for r in records])
    animals, per_animal, grand = animal_average(values, [r["animal"] for r in records])
    _, animal_delta, grand_delta = animal_average(values[:,0]-values[:,1], [r["animal"] for r in records])
    meta = records[0]["meta"]
    t, f = meta["time_s"], meta["frequency_hz"]
    fig, axes = plt.subplots(len(records), 4, squeeze=False,
        figsize=(12, max(3, 2*len(records))), sharex=True, sharey=True, layout="constrained")
    for i, record in enumerate(records):
        for j, (region, route) in enumerate(((0,0), (0,1), (1,0), (1,1))):
            ax = axes[i,j]
            mesh = ax.pcolormesh(t, f, values[i,route,region], shading="auto",
                                 cmap="RdBu_r", vmin=-3, vmax=3, rasterized=True)
            ax.axvline(0, color="black", lw=.6, ls="--")
            if i == 0:
                ax.set_title(f'{("GC","PC")[region]} · {("ortho","retro")[route]}')
            if j == 0:
                ax.set_ylabel(f'{record["animal"]} · {record["session"]}\nFrequency (Hz)', fontsize=6)
            if i == len(records)-1:
                ax.set_xlabel(str(meta["alignment_label"]))
    fig.colorbar(mesh, ax=axes.ravel().tolist(), label="Power (baseline z)", shrink=.5)
    fig.suptitle(f"All-session LFP power · {comparison}")
    save(fig, out, "all_sessions_power")
    fig, axes = plt.subplots(3, 2, squeeze=False, figsize=(7,7), layout="constrained")
    for i, route in enumerate(("Ortho", "Retro", "Ortho − retro")):
        for j, region in enumerate(("GC", "PC")):
            # Difference is formed within session, then averaged within/between animals.
            z = grand[i,j] if i < 2 else grand_delta[j]
            ax = axes[i,j]
            mesh = ax.pcolormesh(t, f, z, shading="auto", cmap="RdBu_r",
                                 vmin=-3, vmax=3, rasterized=True)
            ax.axvline(0, color="black", lw=.6, ls="--")
            ax.set(title=f"{region} · {route}", ylabel="Frequency (Hz)")
            if i == 2:
                ax.set_xlabel(str(meta["alignment_label"]))
    fig.colorbar(mesh, ax=axes.ravel().tolist(), label="Power / difference (z)", shrink=.6)
    fig.suptitle(f"{comparison} · equal-animal mean\n{len(records)} sessions; {len(animals)} animals; descriptive")
    save(fig, out, "animal_mean_power")
    np.savez_compressed(out/"power_summary_values.npz", session_power=values,
        animal_power=per_animal, equal_animal_mean=grand, animal=np.array(animals),
        session=np.array([r["session"] for r in records]),
        session_animal=np.array([r["animal"] for r in records]),
        time_s=t, frequency_hz=f, route=["ortho","retro"], region=["GC","PC"],
        n_animals_per_bin=np.isfinite(per_animal).sum(axis=0),
        animal_difference=animal_delta, equal_animal_difference=grand_delta,
        n_animals_difference=np.isfinite(animal_delta).sum(axis=0))


def plot_phase(records, out, comparison):
    values = np.stack([r["values"] for r in records])
    animals, per_animal, grand = animal_average(values, [r["animal"] for r in records])
    meta = records[0]["meta"]
    t, bands = meta["time_s"], meta["band"]
    fig, axes = plt.subplots(len(records), len(bands), squeeze=False,
        figsize=(3*len(bands), max(3, 1.8*len(records))), sharex=True, sharey=True, layout="constrained")
    for i, record in enumerate(records):
        for j, band in enumerate(bands):
            ax = axes[i,j]
            for k, route in enumerate(("ortho", "retro")):
                ax.plot(t, values[i,k,j], color=COLORS[route], label=route)
            ax.axvline(0, color="black", lw=.6, ls="--")
            ax.set_ylim(0,1)
            if i == 0:
                ax.set_title(str(band))
            if j == 0:
                ax.set_ylabel(f'{record["animal"]} · {record["session"]}\nPhase coherence', fontsize=6)
            if i == len(records)-1:
                ax.set_xlabel(str(meta["alignment_label"]))
    axes[0,0].legend(frameon=False)
    fig.suptitle(f"All-session GC–PC phase coherence · {comparison}")
    save(fig, out, "all_sessions_phase")
    fig, axes = plt.subplots(2, len(bands), squeeze=False,
        figsize=(3*len(bands), 5), sharex=True, layout="constrained")
    for j, band in enumerate(bands):
        for k, route in enumerate(("ortho","retro")):
            axes[0,j].plot(t, per_animal[:,k,j].T, color=COLORS[route], alpha=.25, lw=.6)
            axes[0,j].plot(t, grand[k,j], color=COLORS[route], lw=1.8, label=route)
        delta = per_animal[:,0,j]-per_animal[:,1,j]
        axes[1,j].plot(t, delta.T, color=".6", lw=.6)
        axes[1,j].plot(t, delta.mean(axis=0), color="black", lw=1.8)
        axes[1,j].axhline(0, color=".4", lw=.6)
        axes[0,j].set(title=str(band), ylim=(0,1))
        axes[1,j].set(ylim=(-1,1), xlabel=str(meta["alignment_label"]))
        for ax in axes[:,j]:
            ax.axvline(0, color="black", lw=.6, ls="--")
    axes[0,0].set_ylabel("Phase coherence")
    axes[1,0].set_ylabel("Ortho − retro")
    axes[0,0].legend(frameon=False)
    fig.suptitle(f"{comparison} · {len(records)} sessions; {len(animals)} animals\n"
                 "Thin: animal means; thick: equal-animal mean; descriptive")
    save(fig, out, "animal_mean_phase")
    np.savez_compressed(out/"phase_summary_values.npz", session_coherence=values,
        animal_coherence=per_animal, equal_animal_mean=grand, animal=np.array(animals),
        session=np.array([r["session"] for r in records]),
        session_animal=np.array([r["animal"] for r in records]),
        time_s=t, band=bands, band_limits_hz=meta["band_limits_hz"], route=["ortho","retro"])


def summarize(output, phase_run=None, power_root=None, comparison="route"):
    if comparison not in COMPARISONS:
        raise ValueError("comparison must be route, BA or IA")
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    if phase_run is None and power_root is None:
        raise ValueError("Provide phase_run and/or power_root")
    sources, excluded = [], []
    phase = read_phase(Path(phase_run).resolve(), comparison, sources, excluded) if phase_run is not None else []
    power = read_power(Path(power_root).resolve(), comparison, sources) if power_root is not None else []
    if phase_run is not None and not phase:
        raise ValueError("No completed matched phase conditions")
    if power_root is not None and not power:
        raise ValueError("No session_power.npz found under power root")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".lfp-summary-", dir=output.parent) as tmp:
        staged = Path(tmp)/"summary"; staged.mkdir()
        try:
            with plt.rc_context(STYLE):
                if power:
                    plot_power(power, staged, comparison)
                if phase:
                    plot_phase(phase, staged, comparison)
        finally:
            plt.close("all")
        metadata = dict(comparison=comparison, excluded=excluded,
            weighting="Equal sessions within animal, then equal animals; route = equal BA/IA means",
            uncertainty="Descriptive only; no confidence intervals or significance tests",
            power_sessions=[dict(session=r["session"], animal=r["animal"]) for r in power],
            phase_sessions=[dict(session=r["session"], animal=r["animal"]) for r in phase],
            sources=[dict(path=str(p), sha256=digest(p)) for p in dict.fromkeys(sources)])
        (staged/"summary_manifest.json").write_text(json.dumps(metadata, indent=2))
        if output.exists():
            raise FileExistsError(output)
        staged.rename(output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase-run", type=Path, help="Completed phase run containing status.json")
    parser.add_argument("--power-root", type=Path, help="One run root containing session_power.npz files")
    parser.add_argument("--comparison", choices=tuple(COMPARISONS), default="route")
    parser.add_argument("--output", type=Path, required=True, help="New summary directory")
    args = parser.parse_args()
    print(summarize(args.output, args.phase_run, args.power_root, args.comparison))


if __name__ == "__main__":
    main()
