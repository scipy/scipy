import numpy as np
from scipy.spatial.transform import Rotation, Translation
from numpy.testing import assert_array_almost_equal

def test_rotation_translation_apply():
    # Create a rotation of 90 degrees around the z-axis
    r = Rotation.from_euler('z', 90, degrees=True)
    # Create a translation of (1, 2, 3)
    t = Translation([1, 2, 3])
    
    # Apply rotation and then translation to a point
    point = np.array([1, 0, 0])
    rotated_point = r.apply(point)
    transformed_point = t.apply(rotated_point)
    
    # Expected result after rotation and translation
    expected_point = np.array([1, 3, 3])
    assert_array_almost_equal(transformed_point, expected_point)

def test_translation_rotation_apply():
    # Create a rotation of 90 degrees around the z-axis
    r = Rotation.from_euler('z', 90, degrees=True)
    # Create a translation of (1, 2, 3)
    t = Translation([2, 4, 6])
    
    # Apply translation and then rotation to a point
    point = np.array([1, 0, 0])
    translated_point = t.apply(point)

    # Apply rotation to the translated point (3, 4, 6)
    # Rotating 90 degrees around the z-axis should give (-4, 3, 6)
    transformed_point = r.apply(translated_point)
    
    # Expected result after translation and rotation
    expected_point = np.array([-4, 3, 6])
    assert_array_almost_equal(transformed_point, expected_point)

def test_rotation_and_zero_translation():
    # Create a rotation of 45 degrees around the x-axis
    r1 = Rotation.from_euler('x', 45, degrees=True)
    # Create a translation of (1, 2, 3)
    t1 = Translation([0, 0, 0])
    
    # Apply the first rotation and translation
    point = np.array([1, 1, 1])
    rotated_point = r1.apply(point)
    
    # Apply the second rotation and translation
    translated_rotated_point = r1.apply(t1.apply(point))
    
    # Applying a zero translation should not change the point
    assert_array_almost_equal(rotated_point, translated_rotated_point)

def test_combined_transformations():
    # Create a rotation of 45 degrees around the x-axis
    r1 = Rotation.from_euler('x', 45, degrees=True)
    # Create a translation of (1, 2, 3)
    t1 = Translation([1, 2, 3])
    
    # Create another rotation of 45 degrees around the y-axis
    r2 = Rotation.from_euler('y', 45, degrees=True)
    # Create another translation of (-1, -2, -3)
    t2 = Translation([-1, -2, -3])
    
    # Apply the first rotation and translation
    point = np.array([1, 1, 1])
    transformed_point = t1.apply(r1.apply(point))
    
    # Apply the second rotation and translation
    transformed_point = t2.apply(r2.apply(transformed_point))
    
    # Expected result after combined transformations computed from pydrake
    """
    import numpy as np
    from pydrake.math import RigidTransform, RollPitchYaw
    zero_np = np.zeros(3, dtype=np.float64)

    r1 = RigidTransform(rpy=RollPitchYaw(np.deg2rad(np.array([45.0, 0, 0]))), p=zero_np)
    r2 = RigidTransform(rpy=RollPitchYaw(np.deg2rad(np.array([0, 45.0, 0]))), p=zero_np)
    t1 = RigidTransform(p=np.array([1, 2, 3]))
    t2 = RigidTransform(p=np.array([-1, -2, -3]))

    p = np.array([1.0, 1.0, 1.0])
    p = r1.multiply(p)
    p = t1.multiply(p)
    p = r2.multiply(p)
    p = t2.multiply(p)
    """
    expected_point = np.array([ 3.53553391,  0.        , -1.29289322])
    assert_array_almost_equal(transformed_point, expected_point)

def test_inverse_transformations():
    # Create a rotation of 30 degrees around the z-axis
    r = Rotation.from_euler('z', 30, degrees=True)
    # Create a translation of (1, 2, 3)
    t = Translation([1, 2, 3])
    
    # Apply rotation and translation to a point
    point = np.array([1, 1, 1])
    transformed_point = t.apply(r.apply(point))
    
    # Apply inverse rotation and translation
    inverse_transformed_point = r.inv().apply(t.inv().apply(transformed_point))
    
    # The result should be the original point
    assert_array_almost_equal(inverse_transformed_point, point)
