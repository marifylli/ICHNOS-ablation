"""
CHECK 1 -- is arm B2 still a controlled comparison?

B2 only means anything if the reporter and TIP-TetR modules are identical
across the two repos, so that the only difference is the sensing module. If
the reporter is touched in one repo only, the comparison silently stops being
controlled -- the numbers still come out, they just no longer mean what you
say they mean.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from models import build_snapshots as S


def test_downstream_identical():
    mf = S.verify()
    a, b = mf["downstream_sha"].values()
    assert a == b, (
        "reporter/TIP-TetR DIFFER across repos -- B2 is no longer a controlled "
        "comparison. Do not present the comparison as structural.")


def test_snapshots_unmodified():
    S.verify()   # raises if any SBML was edited after the snapshot


if __name__ == "__main__":
    test_downstream_identical(); test_snapshots_unmodified()
    print("ok: downstream identical; snapshots intact")
