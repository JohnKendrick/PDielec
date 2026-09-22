"""Plot the GaN near-forward modal-pairs validation summary."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyse_results import load_spectra, summarise


ROOT = Path(__file__).resolve().parent
WORKBOOK = ROOT / "Data" / "gan_modal_pairs.xlsx"
FIGURE_DIR = ROOT / "Figures"

CHANNEL_COLOURS = {"e-o": "#0072B2", "o-e": "#D55E00"}
ANGLE_COLOURS = {0.0: "#777777", 2.0: "#009E73"}


def spectrum_label(nac_mode: str, channel: str) -> str:
    """Return a compact legend label."""
    method = "single-q" if nac_mode == "geometry" else "modal pairs"
    initial, final = channel.split("-")
    return rf"{method}, ${initial}\rightarrow {final}$"


def plot_crossed_spectra(ax, frequencies, spectra, angle):
    """Plot geometry and modal-pairs spectra for both crossed channels."""
    mask = (frequencies >= 450.0) & (frequencies <= 780.0)
    for channel in ("e-o", "o-e"):
        colour = CHANNEL_COLOURS[channel]
        for nac_mode, linestyle, linewidth in (
            ("geometry", ":", 2.0),
            ("modal_pairs", "-", 1.7),
        ):
            legend = f"theta_ext={angle:g} {nac_mode} {channel}"
            ax.plot(
                frequencies[mask],
                spectra[legend][mask],
                color=colour,
                linestyle=linestyle,
                linewidth=linewidth,
                label=spectrum_label(nac_mode, channel),
            )

    ax.axvline(485.0, color="#CC79A7", linewidth=1.3, linestyle="--")
    ax.text(
        488.0,
        0.96,
        "experimental finite-$|q|$\nfeature (485 cm$^{-1}$)",
        color="#9C4F7D",
        fontsize=8,
        ha="left",
        va="top",
        transform=ax.get_xaxis_transform(),
    )
    ax.axvline(554.9549, color="0.45", linewidth=1.0, linestyle="-.")
    ax.text(
        558.0,
        0.73,
        r"CRYSTAL23 $\Gamma$ $E_1$(TO)",
        color="0.35",
        fontsize=8,
        rotation=90,
        ha="left",
        va="top",
        transform=ax.get_xaxis_transform(),
    )
    ax.set_xlim(450.0, 780.0)
    ax.set_xlabel(r"Raman shift (cm$^{-1}$)")
    ax.set_title(rf"External incidence {angle:g}$^\circ$")
    ax.grid(axis="y", alpha=0.2)
    if angle == 0.0:
        ax.annotate(
            "modal-pairs crossed\nsignals nearly cancel",
            xy=(650.0, 0.015),
            xytext=(625.0, 1.15),
            arrowprops={"arrowstyle": "->", "color": "0.3", "linewidth": 1.0},
            fontsize=9,
            ha="center",
        )


def plot_integrated_ratios(ax, rows):
    """Plot modal-pairs / geometry integrated-intensity ratios."""
    index = {
        (row["angle_external_deg"], row["nac_mode"], row["channel"]): row
        for row in rows
    }
    channels = ("e-e", "o-o", "e-o", "o-e")
    x = np.arange(len(channels), dtype=float)
    width = 0.36
    for offset, angle in zip((-0.5, 0.5), (0.0, 2.0), strict=True):
        ratios = [
            index[(angle, "modal_pairs", channel)]["integrated_ratio_to_geometry"]
            for channel in channels
        ]
        ax.bar(
            x + offset * width,
            ratios,
            width,
            color=ANGLE_COLOURS[angle],
            label=rf"{angle:g}$^\circ$",
        )

    e_to_o = index[(2.0, "modal_pairs", "e-o")]["integrated_intensity_450_780"]
    o_to_e = index[(2.0, "modal_pairs", "o-e")]["integrated_intensity_450_780"]
    ax.annotate(
        rf"at 2$^\circ$: $I_{{e\to o}}/I_{{o\to e}}={e_to_o / o_to_e:.2f}$",
        xy=(2.5, 0.84),
        xytext=(1.45, 1.08),
        arrowprops={"arrowstyle": "->", "color": "0.3", "linewidth": 1.0},
        fontsize=9,
        ha="left",
    )
    ax.axhline(1.0, color="0.25", linewidth=1.0, linestyle=":")
    tick_labels = []
    for channel in channels:
        initial, final = channel.split("-")
        tick_labels.append(rf"${initial}\rightarrow {final}$")
    ax.set_xticks(x, tick_labels)
    ax.set_ylim(0.0, 1.18)
    ax.set_ylabel("integrated intensity\nmodal pairs / geometry")
    ax.set_title(r"Pair-resolved redistribution, 450–780 cm$^{-1}$")
    ax.text(2.0 - 0.5 * width, 0.025, "≈0", ha="center", va="bottom", fontsize=8)
    ax.text(3.0 - 0.5 * width, 0.025, "≈0", ha="center", va="bottom", fontsize=8)
    ax.legend(frameon=False, title="external angle", ncols=2)
    ax.grid(axis="y", alpha=0.2)


def main() -> None:
    frequencies, spectra = load_spectra(WORKBOOK)
    rows = summarise(frequencies, spectra)

    figure = plt.figure(figsize=(11.0, 7.6), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=(1.0, 0.9))
    axes = [
        figure.add_subplot(grid[0, 0]),
        figure.add_subplot(grid[0, 1]),
        figure.add_subplot(grid[1, :]),
    ]

    plot_crossed_spectra(axes[0], frequencies, spectra, 0.0)
    plot_crossed_spectra(axes[1], frequencies, spectra, 2.0)
    plot_integrated_ratios(axes[2], rows)
    axes[0].set_ylabel("Raman intensity (arb. units)")
    axes[1].set_ylabel("Raman intensity (arb. units)")
    common_ymax = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
    axes[0].set_ylim(0.0, common_ymax)
    axes[1].set_ylim(0.0, common_ymax)

    handles, labels = axes[1].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="outside upper center",
        ncols=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.04),
    )
    figure.suptitle("GaN near-forward Raman: single-q and modal-pairs NAC", fontsize=14)

    FIGURE_DIR.mkdir(exist_ok=True)
    pdf = FIGURE_DIR / "gan_modal_pairs_summary.pdf"
    png = FIGURE_DIR / "gan_modal_pairs_summary.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=220, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote {pdf} and {png}")


if __name__ == "__main__":
    main()
