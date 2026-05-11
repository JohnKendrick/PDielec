PYTEST := pytest

.PHONY: all test verbose fast clean

all: test

test:
	$(PYTEST) .

verbose:
	$(PYTEST) -v .

# Fast tests only (A1-A5, B1, B3, C1, C2, I1-I3 markers)
fast:
	$(PYTEST) test_A_raman_tensor_algebra.py test_B_single_layer.py::TestB1UniformFieldLimit test_B_single_layer.py::TestB3ThicknessScaling test_C_layered_consistency.py::TestC1RamanInactiveLayers test_C_layered_consistency.py::TestC2SubdivisionInvariance test_I_modal_pairs.py::TestI1ModalAmplitudesSumToTotal test_I_modal_pairs.py::TestI2ActivePairSelection test_I_modal_pairs.py::TestI3NACCacheCalledOncePerPair -v

clean:
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null; true
	find . -name "*.pyc" -delete 2>/dev/null; true
