"""Check that qc_table.py takes more input paths than one command line holds.

A run of a few thousand libraries passes qc_table.py more paths than the
kernel's ARG_MAX allows, so the rule writes them to a file. Run with the rule
env: `pixi run --manifest-path workflow/envs/pixi.toml python .tests/unit/test_qc_table_args.py`.
"""

import errno
import os
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).parents[2] / "workflow" / "scripts" / "qc_table.py"


def main():
    """Run qc_table.py on a too-long argument list, then from a file."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        frip = tmp / ("x" * 200) / "e.rep1.frip.tsv"
        frip.parent.mkdir()
        frip.write_text("peak_set\tn_reads\tn_in_peaks\nrep1\t10\t1\n")
        n = os.sysconf("SC_ARG_MAX") // len(str(frip)) + 1
        args = ["--frip", *[str(frip)] * n, "--output", str(tmp / "qc.tsv")]
        try:
            subprocess.run([sys.executable, SCRIPT, *args], check=True)
        except OSError as e:
            assert e.errno == errno.E2BIG, e
        else:
            raise AssertionError(f"{n} paths fit on one command line")
        (tmp / "args.txt").write_text("\n".join(args) + "\n")
        subprocess.run([sys.executable, SCRIPT, f"@{tmp / 'args.txt'}"], check=True)
        assert len((tmp / "qc.tsv").read_text().splitlines()) == n + 1
    print("ok")


if __name__ == "__main__":
    main()
