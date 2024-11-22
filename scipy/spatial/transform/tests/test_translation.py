import numpy as np
from scipy.spatial.transform import Translation

def test_translation_init():
    t = Translation([1, 2, 3])
    assert np.array_equal(t.vector, [1, 2, 3])

def test_translation_add():
    t1 = Translation([1, 2, 3])
    t2 = Translation([4, 5, 6])
    t3 = t1 + t2
    assert np.array_equal(t3.vector, [5, 7, 9])

def test_translation_sub():
    t1 = Translation([4, 5, 6])
    t2 = Translation([1, 2, 3])
    t3 = t1 - t2
    assert np.array_equal(t3.vector, [3, 3, 3])

def test_translation_mul():
    t = Translation([1, 2, 3])
    t2 = t * 2
    assert np.array_equal(t2.vector, [2, 4, 6])

def test_translation_div():
    t = Translation([2, 4, 6])
    t2 = t / 2
    assert np.array_equal(t2.vector, [1, 2, 3])

def test_translation_inv():
    t = Translation([1, 2, 3])
    t_inv = t.inv()
    assert np.array_equal(t_inv.vector, [-1, -2, -3])

def test_translation_magnitude():
    t = Translation([3, 4, 0])
    assert t.magnitude() == 5

def test_translation_approx_equal():
    t1 = Translation([1, 2, 3])
    t2 = Translation([1.00000001, 2.00000001, 3.00000001])
    assert t1.approx_equal(t2)

def test_translation_concatenate():
    t1 = Translation([1, 2, 3])
    t2 = Translation([4, 5, 6])
    t3 = Translation.concatenate([t1, t2])
    assert np.array_equal(t3.vector, [5, 7, 9])

def test_translation_identity():
    t = Translation.identity()
    assert np.array_equal(t.vector, [0, 0, 0])

def test_translation_random():
    t = Translation.random()
    assert t.vector.shape == (3,)
