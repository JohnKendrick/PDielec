"""Finite-q GaN polariton response for the Irmer et al. aperture geometry.

The module supplies both the lossless/damped one-oscillator dielectric model
used in Phys. Rev. B 88, 104303 (2013) and the frequency-dependent dielectric
tensor built by PDielec from the local CRYSTAL23 calculation.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np


LASER_WAVELENGTH_NM = 514.5
SCREEN_DISTANCE_MM = 80.0
Y_VALUES_MM = (2.0, 3.5, 4.5, 6.0, 8.0, 10.0, 12.5, 15.0)
Z_VALUES_MM = (3.0, 4.0, 5.0, 7.0, 10.0)

# Irmer et al., Sec. III B.  Perpendicular means the basal plane and parallel
# means the c axis.  The linewidth is a controlled numerical broadening; the
# paper's analytical dispersion itself was evaluated without damping.
PAPER_EPS_INF_PERP = 5.20
PAPER_EPS_INF_PARALLEL = 5.31
PAPER_TO_E1 = 558.8
PAPER_LO_E1 = 741.0
PAPER_TO_A1 = 531.8
PAPER_LO_A1 = 734.0
DEFAULT_LINEWIDTH_CM1 = 5.0


@dataclass(frozen=True)
class Scan:
    """One experimental aperture and polarization configuration."""

    figure: str
    offset_mm: float
    plane: str
    incident_axis: str
    scattered_axis: str

    @property
    def configuration(self) -> str:
        return f"x({self.incident_axis}{self.scattered_axis})x"


def scans() -> tuple[Scan, ...]:
    """Return the aperture positions shown in Figs. 9--11."""
    result = []
    for y_mm in Y_VALUES_MM:
        result.append(Scan("9", y_mm, "xy", "z", "y"))
    for y_mm in Y_VALUES_MM:
        result.append(Scan("10", y_mm, "xy", "y", "z"))
    for z_mm in Z_VALUES_MM:
        result.append(Scan("11a", z_mm, "xz", "z", "z"))
    for z_mm in Z_VALUES_MM:
        result.append(Scan("11b", z_mm, "xz", "y", "y"))
    return tuple(result)


def paper_permittivity(frequency_cm1: float, linewidth_cm1: float = DEFAULT_LINEWIDTH_CM1):
    """Complex uniaxial dielectric tensor from the paper's TO/LO parameters."""
    v = float(frequency_cm1)

    def component(eps_inf, to, lo):
        denominator = to * to - v * v - 1j * linewidth_cm1 * v
        return eps_inf * (1.0 + (lo * lo - to * to) / denominator)

    eps_perp = component(PAPER_EPS_INF_PERP, PAPER_TO_E1, PAPER_LO_E1)
    eps_parallel = component(PAPER_EPS_INF_PARALLEL, PAPER_TO_A1, PAPER_LO_A1)
    return np.diag([eps_perp, eps_perp, eps_parallel]).astype(complex)


@lru_cache(maxsize=1)
def _dft_permittivity_object():
    """Build the same CRYSTAL23 dielectric object used by PDGui."""
    try:
        from PDielec import Utilities
        from PDielec.HelperRoutines import calculate_dft_permittivity_object
    except ImportError as exc:
        raise RuntimeError(
            "PDielec is required for the dft model; run with the PDielec source "
            "on PYTHONPATH or in an installed PDielec environment"
        ) from exc

    output = Path(__file__).resolve().parents[2] / "Raman" / "raman.log"
    program = Utilities.find_program_from_name(str(output))
    reader = Utilities.get_reader(str(output), program)
    reader.read_output()
    return calculate_dft_permittivity_object(
        reader, sigma=DEFAULT_LINEWIDTH_CM1, eckart=True, mass_definition="Average"
    )


def dft_permittivity(frequency_cm1: float):
    """Complex CRYSTAL23/PDielec dielectric tensor in the crystal frame."""
    return np.asarray(_dft_permittivity_object().calculate(frequency_cm1), dtype=complex)


def epsilon_infinity(model: str) -> np.ndarray:
    """High-frequency principal dielectric constants for optical kinematics."""
    if model == "paper":
        return np.array(
            [PAPER_EPS_INF_PERP, PAPER_EPS_INF_PERP, PAPER_EPS_INF_PARALLEL]
        )
    if model == "dft":
        # Remove the ionic static contribution by reading epsilon infinity from
        # the underlying object; it is the visible-frequency approximation used
        # by PDGui for this calculation.
        obj = _dft_permittivity_object()
        return np.real(np.diag(np.asarray(obj.epsilon_infinity)))
    raise ValueError(f"Unknown dielectric model {model!r}")


def refractive_index(axis: str, model: str) -> float:
    """Visible refractive index for y (ordinary) or z (extraordinary) light."""
    index = 1 if axis == "y" else 2
    return float(np.sqrt(epsilon_infinity(model)[index]))


def aperture_angles(offset_mm: float, scattered_axis: str, model: str):
    """Return external and refracted internal collection angles in degrees."""
    external = np.arctan(float(offset_mm) / SCREEN_DISTANCE_MM)
    n_scattered = refractive_index(scattered_axis, model)
    internal = np.arcsin(np.sin(external) / n_scattered)
    return float(np.degrees(external)), float(np.degrees(internal))


def raman_wavevector(scan: Scan, frequency_cm1: float, model: str) -> np.ndarray:
    """Finite Raman momentum in cycles/cm (the common 2 pi factor is removed)."""
    laser_cm1 = 1.0e7 / LASER_WAVELENGTH_NM
    stokes_cm1 = laser_cm1 - float(frequency_cm1)
    n_laser = refractive_index(scan.incident_axis, model)
    n_stokes = refractive_index(scan.scattered_axis, model)
    external = np.arctan(scan.offset_mm / SCREEN_DISTANCE_MM)
    internal = np.arcsin(np.sin(external) / n_stokes)

    q = np.zeros(3)
    q[0] = n_laser * laser_cm1 - n_stokes * stokes_cm1 * np.cos(internal)
    transverse_axis = 1 if scan.plane == "xy" else 2
    q[transverse_axis] = -n_stokes * stokes_cm1 * np.sin(internal)
    return q


def maxwell_operator(scan: Scan, frequency_cm1: float, model: str) -> np.ndarray:
    """Dimensionless retarded Maxwell operator at self-consistent Raman q."""
    q = raman_wavevector(scan, frequency_cm1, model)
    eps = paper_permittivity(frequency_cm1) if model == "paper" else dft_permittivity(frequency_cm1)
    electromagnetic = (np.dot(q, q) * np.eye(3) - np.outer(q, q)) / frequency_cm1**2
    return electromagnetic - eps


def source_vector(scan: Scan, frequency_cm1: float, model: str) -> np.ndarray:
    """Transverse source selecting the ordinary or extraordinary low branch."""
    q = raman_wavevector(scan, frequency_cm1, model)
    if scan.plane == "xy":
        source = np.array([-q[1], q[0], 0.0])
    else:
        source = np.array([-q[2], 0.0, q[0]])
    return source / np.linalg.norm(source)


def response_point(scan: Scan, frequency_cm1: float, model: str):
    """Return a projected lossy response and a source-independent pole locator."""
    matrix = maxwell_operator(scan, frequency_cm1, model)
    source = source_vector(scan, frequency_cm1, model)
    projected = np.vdot(source, np.linalg.solve(matrix, source))
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    return float(abs(np.imag(projected))), float(1.0 / singular_values[-1])


def frequency_grid() -> np.ndarray:
    """Shared grid resolving both low polariton and parent TO regions."""
    return np.arange(80.0, 590.0 + 0.125, 0.25)


def calculate(model: str):
    """Calculate all spectra and peak summaries for one dielectric model."""
    grid = frequency_grid()
    spectra_rows = []
    summary_rows = []
    for scan in scans():
        response = np.empty_like(grid)
        locator = np.empty_like(grid)
        for index, frequency in enumerate(grid):
            response[index], locator[index] = response_point(scan, frequency, model)

        upper = 557.0 if scan.plane == "xy" else 550.0
        peak_mask = (grid >= 80.0) & (grid <= upper)
        local_indices = np.flatnonzero(peak_mask)
        peak_index = int(local_indices[np.argmax(response[peak_mask])])
        peak_frequency = float(grid[peak_index])
        q = raman_wavevector(scan, peak_frequency, model)
        external_deg, internal_deg = aperture_angles(
            scan.offset_mm, scan.scattered_axis, model
        )
        theta_deg = float(np.degrees(np.arccos(abs(q[2]) / np.linalg.norm(q))))
        summary_rows.append(
            {
                "model": model,
                "figure": scan.figure,
                "configuration": scan.configuration,
                "offset_mm": scan.offset_mm,
                "external_angle_deg": external_deg,
                "internal_collection_angle_deg": internal_deg,
                "polariton_theta_deg": theta_deg,
                "peak_frequency_cm-1": peak_frequency,
                # raman_wavevector omits the common 2 pi so it can be used
                # directly with a Raman shift in cm^-1 in the Maxwell matrix.
                # The paper's dispersion axes report the angular wavevector.
                "q_over_2pi_cm-1": float(np.linalg.norm(q)),
                "q_rad_cm-1": float(2.0 * np.pi * np.linalg.norm(q)),
                "peak_response": float(response[peak_index]),
                "peak_locator": float(locator[peak_index]),
            }
        )
        scale = float(np.max(response))
        for frequency, value, pole in zip(grid, response, locator, strict=True):
            spectra_rows.append(
                {
                    "model": model,
                    "figure": scan.figure,
                    "configuration": scan.configuration,
                    "offset_mm": scan.offset_mm,
                    "frequency_cm-1": float(frequency),
                    "response": float(value),
                    "normalised_response": float(value / scale),
                    "pole_locator": float(pole),
                }
            )
    return spectra_rows, summary_rows


def write_csv(rows: list[dict], path: Path):
    """Write a homogeneous list of dictionaries."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=("paper", "dft", "both"), default="both")
    parser.add_argument("--data-dir", type=Path, default=Path("Data"))
    args = parser.parse_args()

    models = ("paper", "dft") if args.model == "both" else (args.model,)
    all_spectra = []
    all_summary = []
    for model in models:
        spectra_rows, summary_rows = calculate(model)
        all_spectra.extend(spectra_rows)
        all_summary.extend(summary_rows)
    write_csv(all_spectra, args.data_dir / "polariton_response.csv")
    write_csv(all_summary, args.data_dir / "polariton_peaks.csv")
    print(f"Wrote {args.data_dir / 'polariton_response.csv'}")
    print(f"Wrote {args.data_dir / 'polariton_peaks.csv'}")


if __name__ == "__main__":
    main()
