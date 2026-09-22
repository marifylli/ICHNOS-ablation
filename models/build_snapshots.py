"""
Produces IMMUTABLE snapshots of the merged SBML from the two upstream repos.

WHY THIS EXISTS
  If the ablation repo imported the other two repos' code, every upstream
  commit would silently change the results and the comparison would stop
  being reproducible. If instead you copied SBML by hand, one would be
  forgotten and you would be comparing old against new. Here: one script
  produces them, records commit hash + SHA-256, and the runners HALT if
  anything does not match.

  python models/build_snapshots.py
"""
import hashlib, json, subprocess, sys, tempfile, shutil, platform
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
# PIN is the commit the published results were generated from. With
# `--depth 1` and no pin, the MANIFEST recorded a commit that could no longer
# be rebuilt once upstream moved: the hash documented the past but did not
# reproduce it. Set PIN=None only to deliberately move to upstream HEAD, then
# write the new commit back in here.
REPOS = {"FINAL": ("https://github.com/ioanna888/Ichnos-Final-Model.git",
                   "caf14266700071aae4f92d508b6141ff367f04eb"),
         "PULSE": ("https://github.com/ioanna888/Ichnos_PULSE.git",
                   "a4af617f27dc95549e9fd859a0e5794cbc70a5e3")}
VARIANTS = ("ox", "er")
# Files that MUST be identical across repos for arm B2 to be a controlled
# comparison (only the sensing module may differ).
DOWNSTREAM = ("integration/reporter_module_v2.sbml",
              "integration/TIP_TetR_binding.sbml")

BUILD_CODE = """
import sys, io, contextlib, os
sys.path.insert(0, os.getcwd())
from ichnos_core import build_variant_sbml_string
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    s = build_variant_sbml_string(sys.argv[1], save_sbml=False)
sys.stdout.write(s)
"""


def sha(x):
    return hashlib.sha256(x if isinstance(x, bytes) else x.encode()).hexdigest()


def build(work):
    # The solver version decides the numbers, so it belongs in the manifest
    # next to the model hashes -- not only in requirements.txt.
    try:
        import roadrunner, tellurium
        env = {"python": platform.python_version(),
               "libroadrunner": roadrunner.__version__,
               "tellurium": tellurium.__version__}
    except Exception:
        env = {"python": platform.python_version()}
    manifest = {"generated_utc": datetime.now(timezone.utc).isoformat(),
                "environment": env,
                "repos": {}, "snapshots": {}, "downstream_sha": {}}
    helper = work / "_build_one.py"
    helper.write_text(BUILD_CODE)

    for tag, (url, pin) in REPOS.items():
        dst = work / tag
        if pin:
            subprocess.run(["git", "clone", "-q", url, str(dst)], check=True)
            subprocess.run(["git", "-C", str(dst), "checkout", "-q", pin], check=True)
        else:
            subprocess.run(["git", "clone", "-q", "--depth", "1", url, str(dst)], check=True)
        commit = subprocess.run(["git", "-C", str(dst), "rev-parse", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
        if pin and commit != pin:
            raise RuntimeError(f"{tag}: checked out {commit}, expected pin {pin}")
        manifest["repos"][tag] = {"url": url, "commit": commit, "pinned": bool(pin)}
        manifest["downstream_sha"][tag] = {
            f: sha((dst / f).read_bytes()) for f in DOWNSTREAM if (dst / f).exists()}

        for v in VARIANTS:
            r = subprocess.run([sys.executable, str(helper), v],
                               cwd=str(dst / "python"), capture_output=True, text=True)
            if r.returncode != 0 or not r.stdout.strip():
                print(f"  ! {tag}/{v} failed: {r.stderr.strip()[-300:]}")
                continue
            out = HERE / f"merged_{v}_{tag}.sbml"
            # Written and hashed as BYTES. Previously the build hashed
            # r.stdout (a str) while verify() hashed read_text(): on Windows
            # the newline translation made an untouched file look CHANGED.
            out.write_bytes(r.stdout.encode("utf-8"))
            manifest["snapshots"][f"{v}_{tag}"] = {
                "file": out.name, "sha256": sha(out.read_bytes()), "commit": commit}
            print(f"  ok: {out.name}  ({tag} @ {commit[:8]})")

    # The check that makes B2 valid or invalid.
    vals = list(manifest["downstream_sha"].values())
    manifest["downstream_identical"] = (vals[0] == vals[1])
    print("  ok: downstream IDENTICAL (B2 is a controlled comparison)"
          if manifest["downstream_identical"] else
          "  ! WARNING: downstream DIFFERS -- B2 is NOT a controlled comparison")
    (HERE / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def verify():
    """Called by every runner at start-up. Prevents runs against
    out-of-sync or hand-edited models."""
    p = HERE / "MANIFEST.json"
    if not p.exists():
        raise FileNotFoundError("MANIFEST.json missing -- run models/build_snapshots.py")
    mf = json.loads(p.read_text())
    for info in mf["snapshots"].values():
        f = HERE / info["file"]
        if not f.exists():
            raise FileNotFoundError(f"{f.name} missing. Run models/build_snapshots.py")
        if sha(f.read_bytes()) != info["sha256"]:
            raise RuntimeError(f"{info['file']} CHANGED after the snapshot. "
                               "Re-run build_snapshots.py.")
    return mf


def path(variant, tag):
    return HERE / f"merged_{variant}_{tag}.sbml"


if __name__ == "__main__":
    work = Path(tempfile.mkdtemp())
    try:
        build(work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
