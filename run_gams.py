r"""
Runs TriStage_CloudFogEdge.gms for alpha = 0.01, 0.5 and 0.99 from Python (e.g. Spyder)
and then draws the figures with make_figures.py.

GAMS is needed. Either
  (a) the free GAMSPy package, which contains GAMS with a demo licence and CPLEX.
      If GAMS is not found, this script installs it automatically with
      "python -m pip install gamspy" (the demo licence is sufficient for these runs); or
  (b) an installed GAMS system (GAMS Studio). Then set GAMS_DIR below, e.g. r"C:\GAMS\47",
      if it is not found automatically.
"""
import os, sys, glob, shutil, subprocess, runpy

GAMS_DIR = r""                      # optional: folder of an installed GAMS system
ALPHAS = ["0.01", "0.5", "0.99"]
HERE = os.path.dirname(os.path.abspath(__file__))


def find_gams():
    exe = "gams.exe" if os.name == "nt" else "gams"
    cands = [os.path.join(GAMS_DIR, exe)] if GAMS_DIR else []
    try:
        import gamspy_base                      # installed with "pip install gamspy"
        cands.append(os.path.join(gamspy_base.directory, exe))
    except ImportError:
        pass
    w = shutil.which("gams")
    if w:
        cands.append(w)
    for p in [r"C:\GAMS\*\gams.exe", r"C:\GAMS\win64\*\gams.exe", r"C:\Program Files\GAMS*\gams.exe",
              "/opt/gams*/gams", "/opt/gams/*/gams", "/Library/Frameworks/GAMS.framework/Versions/*/Resources/gams"]:
        cands += sorted(glob.glob(p), reverse=True)
    for c in cands:
        if os.path.isfile(c):
            return c
    # no GAMS found: install GAMSPy into the Python that runs this script (works inside Spyder)
    print("GAMS not found - installing the GAMSPy package (GAMS + CPLEX demo) with pip ...", flush=True)
    py = sys.executable
    if os.path.basename(py).lower() == "pythonw.exe":          # Spyder kernel may report pythonw.exe
        py = os.path.join(os.path.dirname(py), "python.exe")
    print("Python used for installation:", py, flush=True)
    ok = False
    for extra in ([], ["--user"]):                             # second try: user site (no admin rights needed)
        r = subprocess.run([py, "-m", "pip", "install", "gamspy"] + extra,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        print((r.stdout or "")[-1500:], (r.stderr or "")[-1500:], flush=True)
        if r.returncode == 0:
            ok = True
            break
    if ok:
        import importlib, site
        importlib.invalidate_caches()
        for p in site.getsitepackages() + [site.getusersitepackages()]:
            if p not in sys.path:
                sys.path.append(p)
        try:
            import gamspy_base
            c = os.path.join(gamspy_base.directory, exe)
            if os.path.isfile(c):
                return c
        except ImportError:
            pass
        sys.exit("GAMSPy was installed. Restart the Spyder kernel (Consoles > Restart kernel) and run this file again.")
    sys.exit("GAMS was not found.\n"
             "Automatic installation of GAMSPy failed (pip messages above). In the Anaconda Prompt (or cmd) type:\n"
             "   \"%s\" -m pip install gamspy\n"
             "or set GAMS_DIR at the top of this file to your GAMS folder." % sys.executable)


def main():
    gams = find_gams()
    print("GAMS:", gams)
    import tempfile
    work = tempfile.mkdtemp(prefix="tristage_")                # ASCII path; GAMS may fail on non-Latin folder names
    shutil.copy(os.path.join(HERE, "TriStage_CloudFogEdge.gms"), work)
    for a in ALPHAS:
        print("solving alpha = %s ..." % a, flush=True)
        r = subprocess.run([gams, "TriStage_CloudFogEdge.gms", "--ALPHA=%s" % a, "lo=2",
                            "o=TriStage_%s.lst" % a, "lf=TriStage_%s.log" % a],
                           cwd=work, capture_output=True, text=True, errors="replace")
        if r.returncode != 0 or not os.path.isfile(os.path.join(work, "results_%s.csv" % a)):
            for f in ("TriStage_%s.log" % a, "TriStage_%s.lst" % a):
                p = os.path.join(work, f)
                if os.path.isfile(p):
                    print("----- end of", f); print(open(p, errors="replace").read()[-2500:])
            sys.exit("GAMS run for alpha = %s failed (return code %s). Files are in %s" % (a, r.returncode, work))
        for f in ("results_%s.csv" % a, "results_%s.gdx" % a, "TriStage_%s.lst" % a):
            shutil.copy(os.path.join(work, f), HERE)
    print("results_<alpha>.csv/.gdx written in", HERE)
    runpy.run_path(os.path.join(HERE, "make_figures.py"), run_name="__main__")


if __name__ == "__main__":
    main()
