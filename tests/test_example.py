import numpy as np
from numpy.testing import assert_allclose

def test_module_imports():
    import castep_adp

def test_bad_array_compare():
    a = np.array([0, 1, 2])
    b = np.array([0, 1, 2])

    assert_allclose(a, b)
