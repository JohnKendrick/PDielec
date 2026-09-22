"""Tests for reader-specific GUI defaults."""

from PDielec.ExperimentOutputReader import ExperimentOutputReader
from PDielec.GenericOutputReader import GenericOutputReader


def test_experimental_reader_disables_optical_symmetrisation_by_default():
    """Experimental tensors should not inherit symmetry from nominal structure data."""
    reader = ExperimentOutputReader([])

    assert reader.default_symmetrise_optical_permittivity is False


def test_generic_reader_enables_optical_symmetrisation_by_default():
    """Electronic-structure readers retain point-group symmetrisation by default."""
    reader = GenericOutputReader([])

    assert reader.default_symmetrise_optical_permittivity is True
