"""Plot the GaN finite-q response in the format of Irmer Figs. 9--11."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyse_verification import load_pdgui_spectra


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
FIGURES = ROOT / "Figures"


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def response_curve(rows, model, figure, configuration, offset):
    selected = [
        row
        for row in rows
        if row["model"] == model
        and row["figure"] == figure
        and row["configuration"] == configuration
        and float(row["offset_mm"]) == float(offset)
    ]
    selected.sort(key=lambda row: float(row["frequency_cm-1"]))
    x = np.asarray([float(row["frequency_cm-1"]) for row in selected])
    y = np.asarray([float(row["normalised_response"]) for row in selected])
    return x, y


def stacked_panel(ax, rows, figure, configuration, offsets, xlim, title):
    """Stack central-ray responses in the style of the experimental figures."""
    spacing = 1.15
    for level, offset in enumerate(offsets):
        baseline = level * spacing
        for model, colour, linestyle, linewidth in (
            ("paper", "#D55E00", "-", 1.6),
            ("dft", "#0072B2", "--", 1.25),
        ):
            x, y = response_curve(rows, model, figure, configuration, offset)
            mask = (x >= xlim[0]) & (x <= xlim[1])
            local = y[mask]
            local = local / np.max(local)
            ax.plot(
                x[mask],
                baseline + local,
                color=colour,
                linestyle=linestyle,
                linewidth=linewidth,
            )
        ax.text(xlim[1] - 4.0, baseline + 0.08, f"{offset:g} mm", fontsize=7.5, ha="right")
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.05, (len(offsets) - 1) * spacing + 1.08)
    ax.set_yticks([])
    ax.set_title(title)
    ax.set_xlabel(r"Raman shift (cm$^{-1}$)")
    ax.grid(axis="x", alpha=0.15)


def plot_spectra(rows):
    figure, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)
    stacked_panel(
        axes[0, 0], rows, "9", "x(zy)x",
        (2, 3.5, 4.5, 6, 8, 10, 12.5, 15), (180, 575),
        r"Fig. 9 geometry: ordinary, $x(zy)x$",
    )
    stacked_panel(
        axes[0, 1], rows, "10", "x(yz)x",
        (2, 3.5, 4.5, 6, 8, 10, 12.5, 15), (80, 575),
        r"Fig. 10 geometry: ordinary, $x(yz)x$",
    )
    stacked_panel(
        axes[1, 0], rows, "11a", "x(zz)x", (3, 4, 5, 7, 10), (250, 570),
        r"Fig. 11(a) geometry: extraordinary, $x(zz)x$",
    )
    stacked_panel(
        axes[1, 1], rows, "11b", "x(yy)x", (3, 4, 5, 7, 10), (250, 570),
        r"Fig. 11(b) geometry: extraordinary, $x(yy)x$",
    )
    handles = [
        plt.Line2D([], [], color="#D55E00", label="published dielectric parameters"),
        plt.Line2D([], [], color="#0072B2", linestyle="--", label="CRYSTAL23/PDielec permittivity"),
    ]
    figure.legend(handles=handles, loc="outside lower center", ncols=2, frameon=False)
    figure.suptitle("GaN central-ray finite-$q$ Maxwell response", fontsize=14)
    FIGURES.mkdir(exist_ok=True)
    for suffix in ("pdf", "png"):
        kwargs = {"dpi": 220} if suffix == "png" else {}
        figure.savefig(FIGURES / f"irmer_response_spectra.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(figure)


def plot_summary(peak_rows, pdgui_frequencies, pdgui_spectra):
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), constrained_layout=True)
    styles = {
        ("9", "x(zy)x"): ("o", r"Fig. 9 $x(zy)x$"),
        ("10", "x(yz)x"): ("s", r"Fig. 10 $x(yz)x$"),
        ("11a", "x(zz)x"): ("^", r"Fig. 11(a) $x(zz)x$"),
        ("11b", "x(yy)x"): ("v", r"Fig. 11(b) $x(yy)x$"),
    }
    for (fig, config), (marker, label) in styles.items():
        for model, colour, linestyle in (
            ("paper", "#D55E00", "-"), ("dft", "#0072B2", "--")
        ):
            selected = [
                row for row in peak_rows
                if row["model"] == model and row["figure"] == fig
                and row["configuration"] == config
            ]
            selected.sort(key=lambda row: float(row["offset_mm"]))
            axes[0].plot(
                [float(row["offset_mm"]) for row in selected],
                [float(row["peak_frequency_cm-1"]) for row in selected],
                marker=marker, color=colour, linestyle=linestyle,
                label=f"{label}, {model}", markersize=4,
            )
    axes[0].set_xlabel("screen displacement $Y$ or $Z$ (mm)")
    axes[0].set_ylabel(r"response maximum (cm$^{-1}$)")
    axes[0].set_title("Finite-$q$ peak moves with aperture")
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=7, frameon=False, ncols=2)

    # Show explicitly that PDGui modal_pairs retains q->0 parent-mode
    # frequencies while the finite-q response moves through the same interval.
    for figure_id, configuration, colour, marker in (
        ("9", "x(zy)x", "#009E73", "o"),
        ("10", "x(yz)x", "#CC79A7", "s"),
    ):
        legends = [
            legend for legend in pdgui_spectra
            if legend.startswith(f"Fig{figure_id} ")
        ]
        offsets, peaks = [], []
        for legend in legends:
            offset = float(legend.split("=")[1].split("mm")[0])
            spectrum = pdgui_spectra[legend]
            mask = (pdgui_frequencies >= 450) & (pdgui_frequencies <= 760)
            peaks.append(float(pdgui_frequencies[mask][np.argmax(spectrum[mask])]))
            offsets.append(offset)
        order = np.argsort(offsets)
        axes[1].plot(
            np.asarray(offsets)[order], np.asarray(peaks)[order],
            marker=marker, color=colour, label=f"PDGui Fig. {figure_id} NAC",
        )
    axes[1].axhline(558.8, color="0.35", linestyle=":", label=r"published $E_1$(TO)")
    axes[1].axhline(741.0, color="0.35", linestyle="--", label=r"published $E_1$(LO)")
    axes[1].set_ylim(430, 765)
    axes[1].set_xlabel("screen displacement $Y$ (mm)")
    axes[1].set_ylabel(r"PDGui spectral maximum (cm$^{-1}$)")
    axes[1].set_title(r"Directional NAC remains a $Gamma$-mode calculation")
    axes[1].grid(alpha=0.2)
    axes[1].legend(fontsize=8, frameon=False)
    for suffix in ("pdf", "png"):
        kwargs = {"dpi": 220} if suffix == "png" else {}
        figure.savefig(FIGURES / f"irmer_verification_summary.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(figure)


def main():
    response_rows = read_rows(DATA / "polariton_response.csv")
    peak_rows = read_rows(DATA / "polariton_peaks.csv")
    frequencies, spectra = load_pdgui_spectra(DATA / "irmer_pdgui.xlsx")
    plot_spectra(response_rows)
    plot_summary(peak_rows, frequencies, spectra)
    print(f"Wrote verification figures under {FIGURES}")


if __name__ == "__main__":
    main()
