"""Summarise and plot the four regenerated GaN powder spectra."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parent


def main():
    """Write band maxima and a plot without normalising individual spectra."""
    with_path = ROOT / "Data" / "gan_powder.xlsx"
    workbook = load_workbook(with_path, read_only=True, data_only=True)
    rows = list(workbook["Powder Raman"].values)
    names = [str(value) for value in rows[0][2:] if value is not None]
    values = np.asarray([row[1:2 + len(names)] for row in rows[1:] if row[1] is not None], dtype=float)
    workbook.close()
    if names != ["none VV", "none VH", "air VV", "air VH"]:
        raise ValueError(f"Unexpected powder scenarios: {names}")
    if not np.all(np.isfinite(values)) or np.any(values[:, 1:] < 0):
        raise ValueError("Powder spectra must be finite and nonnegative")
    frequencies = values[:, 0]
    summaries = []
    for index, name in enumerate(names):
        spectrum = values[:, index + 1]
        for lower, upper in ((100, 200), (500, 590), (600, 780)):
            mask = (frequencies >= lower) & (frequencies <= upper)
            peak = np.argmax(spectrum[mask])
            summaries.append({"scenario": name, "band_cm-1": f"{lower}-{upper}",
                              "maximum_cm-1": frequencies[mask][peak],
                              "peak_intensity": spectrum[mask][peak]})
    output = ROOT / "Data" / "gan_powder_summary.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for ax, matrix in zip(axes, ("none", "air")):
        for index, name in enumerate(names):
            if name.startswith(matrix):
                ax.plot(frequencies, values[:, index + 1], label=name.split()[1])
        ax.set(title=f"Matrix = {matrix}", xlabel=r"Raman shift / cm$^{-1}$",
               ylabel="Powder Raman intensity / a.u.")
        ax.legend()
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    destination = ROOT / "Figures"
    destination.mkdir(exist_ok=True)
    for extension in ("png", "pdf"):
        figure.savefig(destination / f"gan_powder.{extension}", dpi=180)
    plt.close(figure)
    print(f"Wrote {output} and powder figures")


if __name__ == "__main__":
    main()
