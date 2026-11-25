r"""
Regenerates Table 10 (appliance schedules), Table 11 (hourly export/import/storage)
and Figures 6, 8, 9, 10, 11, 12 of the paper from the result files written by
TriStage_CloudFogEdge.gms.

Usage:
  1) gams TriStage_CloudFogEdge.gms --ALPHA=0.01
     gams TriStage_CloudFogEdge.gms --ALPHA=0.5
     gams TriStage_CloudFogEdge.gms --ALPHA=0.99
     (writes results_<alpha>.csv and results_<alpha>.gdx in the same folder)
  2) python make_figures.py   (or press Run in Spyder) - needs only matplotlib;
     it reads the .csv files, so no GAMS tool is needed.
Outputs: fig6.png, fig8.png ... fig12.png, table10_<appliance>.png,
         table11.csv, gams_results.json
"""
import subprocess, csv, io, json, os, shutil, glob, sys, argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ---------------------------------------------------------------- settings (edit here when running in Spyder)
# Folder that contains results_<alpha>.csv / .gdx. Default: the folder of this script.
DATA_DIR = r""
# Only needed if the .csv files are missing and the .gdx files must be read: GAMS system folder,
# e.g. r"C:\GAMS\47". Leave empty to search automatically.
GAMS_DIR = r""

ALPHAS = ["0.01", "0.5", "0.99"]
T = 24
HOURS = list(range(1, T + 1))
COLORS = {"0.01": "#1f77b4", "0.5": "#ff7f0e", "0.99": "#2ca02c"}
APPL = {  # name, allowed clock window [h1,h2) -> periods h1+1..h2
    "CP": ("Computer", (14, 23)), "CO": ("Cooker", (13, 20)),
    "AC": ("Air conditioner", (14, 23)), "WM": ("Washing machine", (9, 17)),
    "PEV": ("Electric vehicle", (16, 24)), "DW": ("Dish washer", (9, 17)),
    "RE": ("Refrigerator", (0, 24)), "VC": ("Vacuum cleaner", (11, 19)),
    "TV": ("Television", (14, 23)), "SD": ("Spin dryer", (9, 17)),
}


def gdxdump():
    """Locate the gdxdump utility that ships with every GAMS installation.
    Search order: --gamsdir argument / GAMSDIR environment variable, PATH,
    the folder of gams(.exe) on PATH, standard install folders, gamspy_base."""
    exe_name = "gdxdump.exe" if os.name == "nt" else "gdxdump"
    cands = []
    if GAMSDIR:
        cands.append(os.path.join(GAMSDIR, exe_name))
    w = shutil.which("gdxdump")
    if w:
        return w
    g = shutil.which("gams")
    if g:
        cands.append(os.path.join(os.path.dirname(os.path.realpath(g)), exe_name))
    pats = [r"C:\GAMS\*\gdxdump.exe", r"C:\GAMS\win64\*\gdxdump.exe",
            r"C:\Program Files\GAMS*\gdxdump.exe", r"C:\Program Files\GAMS\*\gdxdump.exe",
            "/opt/gams/*/gdxdump", "/opt/gams*/gdxdump", "/usr/local/gams*/gdxdump",
            "/Library/Frameworks/GAMS.framework/Versions/*/Resources/gdxdump",
            "/Applications/GAMS*/sysdir/gdxdump", os.path.expanduser("~/gams*/gdxdump")]
    for p in pats:
        cands += sorted(glob.glob(p), reverse=True)   # newest version first
    try:
        import gamspy_base
        cands.append(os.path.join(gamspy_base.directory, exe_name))
    except ImportError:
        pass
    for c in cands:
        if os.path.isfile(c):
            return c
    sys.exit("results_<alpha>.csv not found and gdxdump not found.\n"
             "Re-run the new TriStage_CloudFogEdge.gms (it writes results_<alpha>.csv), or set\n"
             "GAMS_DIR = r\"C:\\GAMS\\47\" at the top of make_figures.py (your GAMS folder).")


GAMSDIR = GAMS_DIR or os.environ.get("GAMSDIR")
_EXE = None


def dump(a, sym):
    global _EXE
    if _EXE is None:
        _EXE = gdxdump()
    f = "results_%s.gdx" % a
    if not os.path.isfile(f):
        sys.exit("results_%s.csv / .gdx not found in %s - run first:  gams TriStage_CloudFogEdge.gms --ALPHA=%s" % (a, os.getcwd(), a))
    out = subprocess.run([_EXE, f, "symb=" + sym, "format=csv"], capture_output=True, text=True).stdout
    return list(csv.reader(io.StringIO(out)))[1:]


def num(s):
    return 0.0 if s in ("Eps", "") else float(s)


def load_csv(a):
    """Read results_<alpha>.csv written by TriStage_CloudFogEdge.gms (no GAMS tool needed)."""
    d = {"demand": [0.0] * T, "on": {}, "ex": [0.0] * T, "im": [0.0] * T, "Ssum": [0.0] * T, "Ysum": [0.0] * T}
    with open("results_%s.csv" % a, newline="") as f:
        for row in csv.DictReader(f):
            sym, t, v = row["symbol"], int(row["t"][1:]) - 1, float(row["value"])
            if sym == "xs":
                d["demand"][t] += v
            elif sym == "xstar":
                if v > 1e-6:
                    d["on"].setdefault(row["a"], {}).setdefault("%s.%s" % (row["r"], row["b"]), []).append(t + 1)
            elif sym in ("ex", "im", "Ssum", "Ysum"):
                d[sym][t] = v
    return d


def load(a):
    if os.path.isfile("results_%s.csv" % a):
        return load_csv(a)
    return load_gdx(a)


def load_gdx(a):
    d = {"demand": [0.0] * T, "on": {}}
    for r, b, t, v in dump(a, "xs"):
        d["demand"][int(t[1:]) - 1] += num(v)
    for r, b, app, t, v in dump(a, "xstar"):
        if num(v) > 1e-6:
            d["on"].setdefault(app, {}).setdefault("%s.%s" % (r, b), []).append(int(t[1:]))
    for sym in ["ex", "im", "Ssum", "Ysum"]:
        v = [0.0] * T
        for row in dump(a, sym):
            v[int(row[0][1:]) - 1] = num(row[-1])
        d[sym] = v
    return d


def bars(ax, key, D, label_fmt="%.1f"):
    w = 0.27
    for k, a in enumerate(ALPHAS):
        xs = [h + (k - 1) * w for h in HOURS]
        bb = ax.bar(xs, D[a][key], w, color=COLORS[a], edgecolor="black", lw=0.4, label="α = %s" % a)
        for rect, val in zip(bb, D[a][key]):
            if val > 1e-6:
                ax.text(rect.get_x() + w / 2, val + 0.25, label_fmt % val, ha="center",
                        va="bottom", rotation=90, fontsize=6.5)
    ax.set_xticks(HOURS)
    ax.set_xlim(0.4, 24.6)
    ax.set_xlabel("Hour")
    ax.grid(axis="y", ls="--", alpha=0.5)
    ax.legend(loc="upper left", fontsize=8)


def find_results():
    """Folder holding results_<alpha>.csv (or .gdx) for all three alphas.
    Looks in DATA_DIR, the folder of this script, the current folder and the
    GAMS Studio workspace (Documents/GAMS/Studio/workspace and sub-folders)."""
    here = os.path.dirname(os.path.abspath(__file__))
    cands = [DATA_DIR, here, os.getcwd()]
    home = os.path.expanduser("~")
    for root in [os.path.join(home, "Documents", "GAMS"), os.path.join(home, "OneDrive", "Documents", "GAMS")]:
        if os.path.isdir(root):
            for dp, dn, fn in os.walk(root):
                if any(f.startswith("results_") for f in fn):
                    cands.append(dp)
    for ext in (".csv", ".gdx"):
        for c in cands:
            if c and all(os.path.isfile(os.path.join(c, "results_%s%s" % (a, ext))) for a in ALPHAS):
                return c
    found = sorted(f for f in os.listdir(here) if f.startswith("results_"))
    sys.exit("The GAMS result files were not found.\n"
             "Needed in the same folder as this script (%s):\n   %s\n"
             "Result files found there: %s\n"
             "Run TriStage_CloudFogEdge.gms (the version in this package) for ALPHA = 0.01, 0.5 and 0.99.\n"
             "GAMS Studio writes the files into its working folder (see the 'Working directory' of the\n"
             "project, usually Documents\\GAMS\\Studio\\workspace); copy them next to this script\n"
             "or set DATA_DIR at the top of this script to that folder."
             % (here, ", ".join("results_%s.csv" % a for a in ALPHAS), ", ".join(found) or "none"))


def main():
    folder = find_results()
    os.chdir(folder)
    print("reading GAMS results from:", folder)
    D = {a: load(a) for a in ALPHAS}
    json.dump(D, open("gams_results.json", "w"), indent=1)
    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 9})

    # Figure 6: hourly consumption
    fig, ax = plt.subplots(figsize=(7.2, 4.5), dpi=200)
    for a, ls, mk in zip(ALPHAS, ["-", "--", ":"], ["o", "s", "^"]):
        ax.plot(HOURS, D[a]["demand"], lw=2, ls=ls, color=COLORS[a], label="α = %s" % a, marker=mk, ms=3.5)
    ax.set_xticks(HOURS); ax.set_xlabel("Hour"); ax.set_ylabel("Energy consumption (kWh)")
    ax.grid(ls="--", alpha=0.5); ax.legend()
    fig.tight_layout(); fig.savefig("fig6.png"); plt.close(fig)

    # Figure 8: total net power supplied to the buildings (sum of y over r,b,t)
    fig, ax = plt.subplots(figsize=(5.0, 3.1), dpi=200)
    tot = [sum(D[a]["Ysum"]) for a in ALPHAS]
    bb = ax.bar(["α = %s" % a for a in ALPHAS], tot, color=[COLORS[a] for a in ALPHAS], edgecolor="black", lw=0.5, width=0.55)
    for r, v in zip(bb, tot):
        ax.text(r.get_x() + r.get_width() / 2, v + 5, "%.2f" % v, ha="center", fontsize=8)
    ax.set_ylabel("Total energy production (kWh)"); ax.set_ylim(0, max(tot) * 1.15)
    ax.grid(axis="y", ls="--", alpha=0.5)
    fig.tight_layout(); fig.savefig("fig8.png"); plt.close(fig)

    # Figure 9: mean hourly export, import and storage level
    fig, ax = plt.subplots(figsize=(4.6, 2.85), dpi=200)
    w = 0.25
    for k, (key, lab, col) in enumerate([("ex", "EX average", "#3a7ca5"), ("im", "IM average", "#f28e2b"),
                                         ("Ssum", "Storage average", "#4e9a50")]):
        vals = [round(sum(D[a][key]) / T + 1e-9, 2) for a in ALPHAS]   # half-up rounding
        bb = ax.bar([i + (k - 1) * w for i in range(3)], vals, w, color=col, edgecolor="black", lw=0.4, label=lab)
        for r, v in zip(bb, vals):
            ax.text(r.get_x() + w / 2, v + 0.3, "%.2f" % v, ha="center", fontsize=6.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(ALPHAS)
    ax.set_xlabel("Consumer comfort level (α)"); ax.set_ylabel("Average energy (kWh)")
    ax.set_ylim(0, 40); ax.grid(axis="y", ls="--", alpha=0.5); ax.legend(fontsize=7, loc="upper left", ncol=3)
    fig.tight_layout(); fig.savefig("fig9.png"); plt.close(fig)

    # Figures 10-12: hourly export, import, storage level
    for name, key, ylab in [("fig10.png", "ex", "Energy exported (kWh)"), ("fig11.png", "im", "Energy imported (kWh)"),
                            ("fig12.png", "Ssum", "Energy stored (kWh)")]:
        fig, ax = plt.subplots(figsize=(10, 4.45), dpi=200)
        bars(ax, key, D)
        ax.set_ylabel(ylab)
        ax.set_ylim(0, max(max(D[a][key]) for a in ALPHAS) * 1.2 + 1)
        fig.tight_layout(); fig.savefig(name); plt.close(fig)

    # Table 11
    with open("table11.csv", "w", newline="") as f:
        w_ = csv.writer(f)
        w_.writerow(["Hour"] + ["%s_a=%s" % (k, a) for a in ALPHAS for k in ("Export", "Import", "Storage")])
        for t in range(T):
            w_.writerow([t + 1] + ["%.2f" % D[a][k][t] for a in ALPHAS for k in ("ex", "im", "Ssum")])

    # Table 10: one schedule grid per appliance (rows: region, building, alpha; columns: hours)
    blds = sorted({rb for a in ALPHAS for app in D[a]["on"] for rb in D[a]["on"][app]})
    for app, (title, (h1, h2)) in APPL.items():
        win = set(range(h1 + 1, h2 + 1))
        nrow = len(blds) * len(ALPHAS)
        fig, ax = plt.subplots(figsize=(3.45, 2.5), dpi=300)
        for i, rb in enumerate(blds):
            for j, a in enumerate(ALPHAS):
                y = nrow - 1 - (i * len(ALPHAS) + j)
                on = set(D[a]["on"].get(app, {}).get(rb, []))
                for t in HOURS:
                    ax.add_patch(Rectangle((t - 1, y), 1, 1, facecolor="#808080" if t in on else "white",
                                           edgecolor="black" if t in win else "#bbbbbb", lw=0.5 if t in win else 0.25))
                ax.text(-0.3, y + 0.5, a, ha="right", va="center", fontsize=3.2)
            ax.text(-3.2, nrow - (i + 0.5) * len(ALPHAS), rb.replace("r", "R").replace(".b", ", B"),
                    ha="right", va="center", fontsize=3.6)
            ax.axhline(nrow - (i + 1) * len(ALPHAS), color="black", lw=0.7)
        for t in HOURS:
            ax.text(t - 0.5, nrow + 0.4, str(t), ha="center", va="bottom", fontsize=3.4)
        ax.set_xlim(-6.5, 24.1); ax.set_ylim(-1.8, nrow + 1.6); ax.axis("off")
        ax.text(8.8, -1.0, title, ha="center", va="center", fontsize=6)
        fig.tight_layout(pad=0.1); fig.savefig("table10_%s.png" % app); plt.close(fig)

    for a in ALPHAS:
        d = D[a]
        pk = max(d["demand"])
        print("alpha=%s peak=%.2f at t=%d  export=%.2f  import=%.2f  sum(y)=%.2f  mean storage=%.2f"
              % (a, pk, d["demand"].index(pk) + 1, sum(d["ex"]), sum(d["im"]), sum(d["Ysum"]), round(sum(d["Ssum"]) / T + 1e-9, 2)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamsdir", help="GAMS system folder that contains gdxdump")
    args, _ = ap.parse_known_args()
    if args.gamsdir:
        GAMSDIR = args.gamsdir
    main()
