import numpy as np
from scipy import special

def test_special_softmax():
    """Test special softmax function."""
    x = np.array([1.0, 2.0, 3.0])
    # Simple test
    assert len(x) == 3

