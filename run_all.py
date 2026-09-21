"""
Run the complete ablation study end to end, in the correct order.

  python run_all.py              # everything
  python run_all.py --no-build   # reuse the existing model snapshots

Order matters: snapshots -> checks -> arms. If a check fails, the arms do not
run, because their results would not mean what the wiki says they mean.
Arm B1 is expected to stop cleanly until the Delaunay data file is added.
"""
import subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent

BUILD  = [["models/build_snapshots.py"]]
CHECKS = [["checks/test_downstream_identical.py"],
          ["checks/test_protocol_invariance.py"]]
ARMS   = [["arms/run_A_lesions.py", "PULSE", "ox"],
          ["arms/run_A_lesions.py", "PULSE", "er"],
          ["arms/run_A3_sensitivity.py", "PULSE", "ox"],
          ["arms/run_A3_sensitivity.py", "PULSE", "er"],
          ["arms/run_A4_kd_spec.py", "PULSE", "ox"],
          ["arms/run_B2_structure.py", "ox"],
          ["arms/run_B2_structure.py", "er"]]
OPTIONAL = [["arms/run_B1_aicc.py"]]   # halts until data/ is populated


def step(cmd, required=True):
    t0 = time.time()
    print(f"\n{'=' * 72}\n>> {' '.join(cmd)}\n{'=' * 72}", flush=True)
    r = subprocess.run([sys.executable, *cmd], cwd=ROOT)
    dt = time.time() - t0
    if r.returncode != 0:
        if required:
            print(f"\n!! FAILED ({dt:.0f}s): {' '.join(cmd)} -- stopping.")
            sys.exit(r.returncode)
        print(f"\n-- skipped ({dt:.0f}s): {' '.join(cmd)} (expected until its data exists)")
        return False
    print(f"-- ok ({dt:.0f}s)")
    return True


if __name__ == "__main__":
    if "--no-build" not in sys.argv:
        for c in BUILD:
            step(c)
    for c in CHECKS:
        step(c)
    for c in ARMS:
        step(c)
    for c in OPTIONAL:
        step(c, required=False)
    print(f"\nAll done. Results in {ROOT / 'results'}")
