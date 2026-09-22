"""Tests D1-D4: System class layer management and geometry.

Plan sections D1 (add/del/get layer), D2 (substrate/superstrate),
D3 (get_layers_boundaries), D4 (reversed_system).
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

import numpy as np
import pytest
from PDielec.GTMcore import Layer, System, TransferMatrixSystem


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _iso_eps(n):
    """Return a constant isotropic permittivity function for refractive index n."""
    eps = complex(n) ** 2
    return lambda freq_cm1: eps * np.eye(3, dtype=complex)


def _make_layer(thickness_m, n=1.0):
    """Create a Layer with coherent mode flags."""
    layer = Layer(thickness=thickness_m, epsilon=_iso_eps(n))
    layer.inCoherentIntensity = False
    layer.inCoherentPhase = False
    layer.inCoherentAveragePhase = False
    layer.inCoherentThick = False
    layer.SMatrix = None
    return layer


# ---------------------------------------------------------------------------
# D1: add_layer / del_layer / get_layer
# ---------------------------------------------------------------------------

class TestD1LayerManagement:
    """D1: Adding, removing, and retrieving layers from a System."""

    def test_empty_system_has_no_layers(self):
        sys_ = System()
        assert sys_.get_all_layers() == []

    def test_add_one_layer(self):
        sys_ = System()
        L = _make_layer(1e-6)
        sys_.add_layer(L)
        assert len(sys_.get_all_layers()) == 1

    def test_get_layer_returns_correct_object(self):
        sys_ = System()
        L = _make_layer(1e-6, n=1.5)
        sys_.add_layer(L)
        assert sys_.get_layer(0) is L

    def test_del_layer_removes_it(self):
        sys_ = System()
        L = _make_layer(1e-6)
        sys_.add_layer(L)
        sys_.del_layer(0)
        assert sys_.get_all_layers() == []

    def test_del_invalid_position_does_not_crash(self):
        sys_ = System()
        # Should log an error but not raise
        sys_.del_layer(5)

    def test_add_two_layers(self):
        sys_ = System()
        L1 = _make_layer(1e-6)
        L2 = _make_layer(2e-6)
        sys_.add_layer(L1)
        sys_.add_layer(L2)
        assert len(sys_.get_all_layers()) == 2

    def test_layer_order_preserved(self):
        sys_ = System()
        L1 = _make_layer(1e-6)
        L2 = _make_layer(2e-6)
        sys_.add_layer(L1)
        sys_.add_layer(L2)
        assert sys_.get_layer(0) is L1
        assert sys_.get_layer(1) is L2


# ---------------------------------------------------------------------------
# D2: get_substrate / get_superstrate
# ---------------------------------------------------------------------------

class TestD2SubstrateSuperstate:
    """D2: Setting and retrieving substrate and superstrate."""

    def test_get_substrate_returns_set_substrate(self):
        sub = _make_layer(1e-3, n=1.5)
        sys_ = System(substrate=sub)
        assert sys_.get_substrate() is sub

    def test_get_superstrate_returns_set_superstrate(self):
        sup = _make_layer(1e-3, n=1.0)
        sys_ = System(superstrate=sup)
        assert sys_.get_superstrate() is sup

    def test_set_substrate_round_trip(self):
        sys_ = System()
        sub = _make_layer(1e-3, n=2.0)
        sys_.set_substrate(sub)
        assert sys_.get_substrate() is sub

    def test_set_superstrate_round_trip(self):
        sys_ = System()
        sup = _make_layer(1e-3, n=1.0)
        sys_.set_superstrate(sup)
        assert sys_.get_superstrate() is sup


# ---------------------------------------------------------------------------
# D3: get_layers_boundaries
# ---------------------------------------------------------------------------

class TestD3LayerBoundaries:
    """D3: Boundary positions are accumulated correctly."""

    def test_no_layers_boundaries(self):
        sup = _make_layer(1e-4, n=1.0)
        sub = _make_layer(2e-4, n=1.5)
        sys_ = System(substrate=sub, superstrate=sup)
        zn = sys_.get_layers_boundaries()
        # [superstrate bottom, 0, substrate top] = 3 elements minimum
        assert len(zn) >= 3

    def test_one_layer_boundary_positions(self):
        sup = _make_layer(1e-4, n=1.0)
        sub = _make_layer(2e-4, n=1.5)
        d = 5e-7
        L = _make_layer(d, n=2.0)
        sys_ = System(substrate=sub, superstrate=sup, layers=[L])
        zn = sys_.get_layers_boundaries()
        # zn[1] = 0 (entrance surface), zn[2] = d (exit of finite layer)
        assert abs(zn[1]) < 1e-15
        assert abs(zn[2] - d) < 1e-15

    def test_two_layers_cumulative_boundaries(self):
        sup = _make_layer(1e-4)
        sub = _make_layer(1e-4)
        d1, d2 = 1e-7, 3e-7
        sys_ = System(substrate=sub, superstrate=sup,
                      layers=[_make_layer(d1), _make_layer(d2)])
        zn = sys_.get_layers_boundaries()
        assert abs(zn[2] - d1) < 1e-15
        assert abs(zn[3] - (d1 + d2)) < 1e-15


# ---------------------------------------------------------------------------
# D4: reversed_system
# ---------------------------------------------------------------------------

class TestD4ReversedSystem:
    """D4: reversed_system swaps substrate/superstrate and reverses layer order."""

    def test_reversed_same_number_of_layers(self):
        sup = _make_layer(1e-3, n=1.0)
        sub = _make_layer(1e-3, n=1.5)
        L1 = _make_layer(1e-7, n=2.0)
        L2 = _make_layer(2e-7, n=3.0)
        sys_ = TransferMatrixSystem(substrate=sub, superstrate=sup, layers=[L1, L2])
        rev = sys_.reversed_system()
        assert len(rev.get_all_layers()) == 2

    def test_reversed_substrate_becomes_superstrate(self):
        sup = _make_layer(1e-3, n=1.0)
        sub = _make_layer(1e-3, n=1.5)
        sys_ = TransferMatrixSystem(substrate=sub, superstrate=sup)
        rev = sys_.reversed_system()
        assert rev.get_superstrate() is sub
        assert rev.get_substrate() is sup

    def test_layer_order_is_reversed(self):
        sup = _make_layer(1e-3)
        sub = _make_layer(1e-3)
        L1 = _make_layer(1e-7, n=2.0)
        L2 = _make_layer(2e-7, n=3.0)
        sys_ = TransferMatrixSystem(substrate=sub, superstrate=sup, layers=[L1, L2])
        rev = sys_.reversed_system()
        layers = rev.get_all_layers()
        assert layers[0] is L2
        assert layers[1] is L1

    def test_reversed_is_same_type(self):
        sup = _make_layer(1e-3)
        sub = _make_layer(1e-3)
        sys_ = TransferMatrixSystem(substrate=sub, superstrate=sup)
        rev = sys_.reversed_system()
        assert isinstance(rev, TransferMatrixSystem)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
