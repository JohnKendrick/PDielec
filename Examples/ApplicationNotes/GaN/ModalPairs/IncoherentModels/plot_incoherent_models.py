"""Plot the GaN incoherent final-state model validation."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyse_incoherent_models import load_spectra, summarise


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "Figures"
CHANNEL_COLOURS = {"e-o": "#0072B2", "o-e": "#D55E00"}


def legend(angle, model, tolerance, channel):
    """Return the workbook legend for one scenario."""
    return (
        f"theta_ext={angle:g} model={model} tol={tolerance:g} "
        f"channel={channel}"
    )


def channel_tex(channel):
    """Return a mathtext channel label."""
    initial, final = channel.split("-")
    return rf"${initial}\rightarrow {final}$"


def plot_spectra(ax, frequencies, spectra, angle):
    """Compare bulk and local spectra for the crossed channels."""
    mask = (frequencies >= 450.0) & (frequencies <= 780.0)
    for channel in ("e-o", "o-e"):
        colour = CHANNEL_COLOURS[channel]
        ax.plot(
            frequencies[mask],
            spectra[legend(angle, "local", 90.0, channel)][mask],
            color=colour,
            linestyle=":",
            linewidth=2.2,
            label=f"local, {channel_tex(channel)}",
        )
        ax.plot(
            frequencies[mask],
            spectra[legend(angle, "bulk", 90.0, channel)][mask],
            color=colour,
            linestyle="-",
            linewidth=1.8,
            label=f"bulk 90°, {channel_tex(channel)}",
        )
        if angle == 2.0:
            ax.plot(
                frequencies[mask],
                spectra[legend(angle, "bulk", 0.0, channel)][mask],
                color=colour,
                linestyle="--",
                linewidth=1.2,
                label=f"bulk 0°, {channel_tex(channel)}",
            )
    ax.set_xlim(450.0, 780.0)
    ax.set_xlabel(r"Raman shift (cm$^{-1}$)")
    ax.set_ylabel("Raman intensity (arb. units)")
    ax.set_title(rf"External incidence {angle:g}$^\circ$: crossed channels")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(
        frameon=False,
        fontsize=8,
        ncols=2,
        loc="upper right" if angle == 0.0 else "upper left",
    )


def indexed(rows):
    """Index summary rows by scenario parameters."""
    return {
        (
            row["angle_external_deg"],
            row["model"],
            row["q_tolerance_deg"],
            row["channel"],
        ): row
        for row in rows
    }


def intensity(index, angle, model, tolerance, channel):
    """Return one integrated polar-band intensity."""
    return index[(angle, model, tolerance, channel)][
        "integrated_intensity_450_780"
    ]


def plot_model_ratio(ax, rows):
    """Plot bulk/local ratios at the 90 degree acceptance setting."""
    index = indexed(rows)
    channels = ("e-e", "o-o", "e-o", "o-e")
    x = np.arange(4, dtype=float)
    width = 0.36
    for offset, angle, colour in (
        (-0.5, 0.0, "#777777"),
        (0.5, 2.0, "#009E73"),
    ):
        values = [
            intensity(index, angle, "bulk", 90.0, channel)
            / intensity(index, angle, "local", 90.0, channel)
            for channel in channels
        ]
        ax.bar(
            x + offset * width,
            values,
            width,
            color=colour,
            label=rf"{angle:g}$^\circ$",
        )
    ax.axhline(1.0, color="0.25", linestyle=":", linewidth=1.0)
    ax.set_xticks(x, [channel_tex(channel) for channel in channels])
    ax.set_ylim(0.0, 1.08)
    ax.set_ylabel("bulk / local integrated intensity")
    ax.set_title("Effect of final-state model at 90° tolerance")
    ax.legend(frameon=False, title="external angle")
    ax.grid(axis="y", alpha=0.2)


def plot_tolerance_ratio(ax, rows):
    """Show that tolerance affects bulk but not local incoherent spectra."""
    index = indexed(rows)
    channels = ("e-e", "o-o", "e-o", "o-e")
    x = np.arange(4, dtype=float)
    width = 0.36
    bulk = [
        intensity(index, 2.0, "bulk", 0.0, channel)
        / intensity(index, 2.0, "bulk", 90.0, channel)
        for channel in channels
    ]
    local = [
        intensity(index, 2.0, "local", 0.0, channel)
        / intensity(index, 2.0, "local", 90.0, channel)
        for channel in channels
    ]
    ax.bar(x - 0.5 * width, bulk, width, color="#CC79A7", label="bulk")
    ax.bar(x + 0.5 * width, local, width, color="#56B4E9", label="local")
    ax.axhline(1.0, color="0.25", linestyle=":", linewidth=1.0)
    ax.set_yscale("log")
    ax.set_ylim(1.0e-12, 2.0)
    ax.set_xticks(x, [channel_tex(channel) for channel in channels])
    ax.set_ylabel("integrated intensity, tol 0° / tol 90°")
    ax.set_title("Tolerance control at 2° external incidence")
    ax.legend(frameon=False)
    ax.grid(axis="y", which="both", alpha=0.2)


def main():
    frequencies, spectra = load_spectra()
    rows = summarise(frequencies, spectra)
    figure, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), constrained_layout=True)

    plot_spectra(axes[0, 0], frequencies, spectra, 0.0)
    plot_spectra(axes[0, 1], frequencies, spectra, 2.0)
    plot_model_ratio(axes[1, 0], rows)
    plot_tolerance_ratio(axes[1, 1], rows)

    common_ymax = max(axes[0, 0].get_ylim()[1], axes[0, 1].get_ylim()[1])
    axes[0, 0].set_ylim(0.0, common_ymax)
    axes[0, 1].set_ylim(0.0, common_ymax)
    figure.suptitle("GaN modal-pairs: incoherent final-state models", fontsize=14)

    FIGURE_DIR.mkdir(exist_ok=True)
    pdf = FIGURE_DIR / "gan_incoherent_models.pdf"
    png = FIGURE_DIR / "gan_incoherent_models.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=220, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote {pdf} and {png}")


if __name__ == "__main__":
    main()
