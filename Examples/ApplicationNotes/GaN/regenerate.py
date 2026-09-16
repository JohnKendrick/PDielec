"""Regenerate the GaN application-note workbooks, summaries and figures.

Use the Python environment containing PDielec. All paths are relative to this
file, so the driver can be launched from any directory. Existing generated
products are overwritten; the CRYSTAL calculation is not rerun.
"""

import argparse
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
JOBS = {
    "powder": ("Powder", "gan_powder", ["analyse_powder.py"]),
    "modal": ("ModalPairs", "gan_modal_pairs", ["analyse_results.py", "plot_results.py"]),
    "incoherent": ("ModalPairs/IncoherentModels", "gan_incoherent_models",
                   ["analyse_incoherent_models.py", "plot_incoherent_models.py"]),
    "irmer": ("ModalPairs/Irmer2013", "irmer_pdgui",
              ["polariton_response.py", "analyse_verification.py", "plot_verification.py"]),
}


def main():
    """Run selected examples, stopping immediately if any subprocess fails."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("examples", nargs="*", help="powder, modal, incoherent, irmer (default: all)")
    parser.add_argument("--cpus", type=int, default=2)
    args = parser.parse_args()
    selected = args.examples or list(JOBS)
    if args.cpus < 1 or any(name not in JOBS for name in selected):
        parser.error("choose known examples and a positive --cpus value")
    env = os.environ.copy()
    checkout = ROOT.parents[2]
    if (checkout / "PDielec" / "pdgui.py").is_file():
        env["PYTHONPATH"] = str(checkout) + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("QT_QPA_PLATFORM", "offscreen")
    env.setdefault("MPLBACKEND", "Agg")
    for name in selected:
        directory, stem, analyses = JOBS[name]
        cwd = ROOT / directory
        (cwd / "Data").mkdir(exist_ok=True)
        workbook = cwd / "Data" / f"{stem}.xlsx"
        commands = [[sys.executable, "-m", "PDielec.pdgui", "-cpus", str(args.cpus),
                     "-nosplash", "-exit", "-script", f"{stem}.py",
                     "-spreadsheet", str(workbook)]]
        commands.extend([sys.executable, script] for script in analyses)
        for command in commands:
            print(f"[{name}] {' '.join(command)}", flush=True)
            subprocess.run(command, cwd=cwd, env=env, check=True)


if __name__ == "__main__":
    main()
