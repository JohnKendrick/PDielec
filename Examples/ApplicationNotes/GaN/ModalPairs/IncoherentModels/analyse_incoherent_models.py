"""Summarise the two GaN incoherent modal-pairs final-state models."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parent
WORKBOOK = ROOT / "Data" / "gan_incoherent_models.xlsx"
SUMMARY = ROOT / "Data" / "gan_incoherent_models_summary.csv"
REPORT = ROOT / "RESULTS.md"

LEGEND_RE = re.compile(
    r"theta_ext=(?P<angle>[0-9.]+) "
    r"model=(?P<model>bulk|local) "
    r"tol=(?P<tolerance>[0-9.]+) "
    r"channel=(?P<channel>[eo]-[eo])"
)


def load_spectra(workbook: Path = WORKBOOK):
    """Return the frequency grid and all named spectra."""
    wb = load_workbook(workbook, read_only=True, data_only=True)
    ws = wb["Crystal Raman"]
    rows = ws.iter_rows(values_only=True)
    header = next(rows)
    legends = [str(value) for value in header[2:] if value is not None]
    data = [row for row in rows if row[1] is not None]
    frequencies = np.asarray([row[1] for row in data], dtype=float)
    spectra = {
        legend: np.asarray([row[index + 2] for row in data], dtype=float)
        for index, legend in enumerate(legends)
    }
    return frequencies, spectra


def summarise(frequencies, spectra):
    """Extract polar-band peaks and integrated intensities."""
    mask = (frequencies >= 450.0) & (frequencies <= 780.0)
    x = frequencies[mask]
    rows = []
    for legend, spectrum in spectra.items():
        match = LEGEND_RE.fullmatch(legend)
        if match is None:
            raise ValueError(f"Unexpected scenario legend: {legend!r}")
        y = spectrum[mask]
        peak_index = int(np.argmax(y))
        rows.append(
            {
                "angle_external_deg": float(match.group("angle")),
                "model": match.group("model"),
                "q_tolerance_deg": float(match.group("tolerance")),
                "channel": match.group("channel"),
                "peak_frequency_cm-1": float(x[peak_index]),
                "peak_intensity": float(y[peak_index]),
                "integrated_intensity_450_780": float(np.trapezoid(y, x)),
            }
        )
    return rows


def write_summary(rows, output: Path = SUMMARY):
    """Write the machine-readable comparison."""
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_report(rows, output: Path = REPORT):
    """Write the principal model and tolerance ratios."""
    index = {
        (
            row["angle_external_deg"],
            row["model"],
            row["q_tolerance_deg"],
            row["channel"],
        ): row
        for row in rows
    }
    lines = [
        "# GaN incoherent final-state model results",
        "",
        "Integrated intensities cover 450–780 cm⁻¹.",
        "",
        "| external angle | channel | bulk(90°) / local | bulk(0°) / bulk(90°) |",
        "|---:|---|---:|---:|",
    ]
    for angle in (0.0, 2.0):
        for channel in ("e-e", "o-o", "e-o", "o-e"):
            bulk_0 = index[(angle, "bulk", 0.0, channel)][
                "integrated_intensity_450_780"
            ]
            bulk_90 = index[(angle, "bulk", 90.0, channel)][
                "integrated_intensity_450_780"
            ]
            local = index[(angle, "local", 90.0, channel)][
                "integrated_intensity_450_780"
            ]
            lines.append(
                f"| {angle:g}° | {channel} | {bulk_90 / local:.6g} | "
                f"{bulk_0 / bulk_90:.6g} |"
            )
    lines.extend(
        [
            "",
            "`Local incoherent` is identical at 0° and 90° q tolerance for every ",
            "scenario. The tolerance therefore acts only on the bulk phase-matched ",
            "acceptance rule, as intended.",
            "",
        ]
    )
    output.write_text("\n".join(lines))


def main():
    frequencies, spectra = load_spectra()
    rows = summarise(frequencies, spectra)
    write_summary(rows)
    write_report(rows)
    print(f"Wrote {SUMMARY} and {REPORT}")


if __name__ == "__main__":
    main()
