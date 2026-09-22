"""Tests for Crystal Raman Porto notation parsing."""

import pytest
import numpy as np

from PDielec.GUI.CrystalScenarioTab import _parse_porto_notation


@pytest.mark.parametrize(
    ("porto", "hkl", "incident", "detected", "side"),
    [
        ("x(yy)x", [1, 0, 0], "p", "p", "superstrate"),
        ("-x(yy)x", [1, 0, 0], "p", "p", "superstrate"),
        ("x(yz)-x", [1, 0, 0], "p", "s", "superstrate"),
        (r"x(yz)\bar{x}", [1, 0, 0], "p", "s", "superstrate"),
        ("x(yz)+x", [1, 0, 0], "p", "s", "substrate"),
        ("-x(yz)-x", [1, 0, 0], "p", "s", "substrate"),
        ("z(xx)z", [0, 0, 1], "p", "p", "superstrate"),
        ("z(xy)barz", [0, 0, 1], "p", "s", "superstrate"),
    ],
)
def test_parse_porto_notation(porto, hkl, incident, detected, side):
    """Principal-axis Porto notation maps to current Crystal Raman settings."""
    mapped = _parse_porto_notation(porto)

    assert mapped["hkl"] == hkl
    assert mapped["incident_pol"] == incident
    assert mapped["detected_pol"] == detected
    assert mapped["collection_side"] == side
    assert mapped["angle_of_incidence"] == 0.0
    assert mapped["global_azimuthal_angle"] == 0.0


def test_parse_porto_notation_compound_polarisation():
    """Compound Porto polarisations map to normalised p/s Jones vectors."""
    mapped = _parse_porto_notation("z(x+y,x+y)z")

    assert mapped["hkl"] == [0, 0, 1]
    np.testing.assert_allclose(mapped["incident_pol"], np.array([1.0, 1.0]) / np.sqrt(2.0))
    np.testing.assert_allclose(mapped["detected_pol"], np.array([1.0, 1.0]) / np.sqrt(2.0))
    assert mapped["collection_side"] == "superstrate"


@pytest.mark.parametrize("porto", ["x(xx)x", "x(yy)z", "z(x+y)z", "bad", ""])
def test_parse_porto_notation_rejects_unsupported_forms(porto):
    """Invalid and non-collinear forms are rejected explicitly."""
    with pytest.raises(ValueError):
        _parse_porto_notation(porto)
