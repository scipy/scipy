# cython: cpow=True

import numpy as np
from ._rotation import Rotation, compose_quat

cimport numpy as np
cimport cython

np.import_array()


cdef _create_skew_matrix(vec):
    """Create skew-symmetric (aka cross-product) matrix for stack of vectors."""
    result = np.zeros((len(vec), 3, 3))
    result[:, 0, 1] = -vec[:, 2]
    result[:, 0, 2] = vec[:, 1]
    result[:, 1, 0] = vec[:, 2]
    result[:, 1, 2] = -vec[:, 0]
    result[:, 2, 0] = -vec[:, 1]
    result[:, 2, 1] = vec[:, 0]
    return result


cdef _compute_se3_exp_translation_transform(rot_vec):
    """Compute the transformation matrix from the se3 translation part to SE3
    translation.

    The transformation matrix depends on the rotation vector.
    """
    angle = np.linalg.norm(rot_vec, axis=1)
    mask = angle < 1e-3

    k1 = np.empty(len(rot_vec))
    k1[mask] = 0.5 - angle[mask]**2 / 24 + angle[mask]**4 / 720
    k1[~mask] = (1.0 - np.cos(angle[~mask])) / angle[~mask]**2

    k2 = np.empty(len(rot_vec))
    k2[mask] = 1/6 - angle[mask]**2 / 120 + angle[mask]**4 / 5040
    k2[~mask] = (angle[~mask] - np.sin(angle[~mask])) / angle[~mask]**3

    s = _create_skew_matrix(rot_vec)

    return np.identity(3) + k1[:, None, None] * s + k2[:, None, None] * s @ s


cdef _compute_se3_log_translation_transform(rot_vec):
    """Compute the transformation matrix from SE3 translation to the se3
    translation part.

    It is the inverse of `_compute_se3_exp_translation_transform` in a closed
    analytical form.
    """
    angle = np.linalg.norm(rot_vec, axis=1)
    mask = angle < 1e-3

    k = np.empty(len(rot_vec))
    k[mask] = 1/12 + angle[mask]**2 / 720 + angle[mask]**4 / 30240
    k[~mask] = (1 - 0.5 * angle[~mask] / np.tan(0.5 * angle[~mask])) / angle[~mask]**2

    s = _create_skew_matrix(rot_vec)

    return np.identity(3) - 0.5 * s + k[:, None, None] * s @ s


def normalize_dual_quaternion(dual_quat):
    """Normalize dual quaternion."""
    real, dual = _normalize_dual_quaternion(
        dual_quat[..., :4], dual_quat[..., 4:])
    return np.concatenate((real, dual), axis=-1)


cdef _normalize_dual_quaternion(real_part, dual_part):
    """Ensure that unit norm of the dual quaternion.

    The norm is a dual number and must be 1 + 0 * epsilon, which means that
    the real quaternion must have unit norm and the dual quaternion must be
    orthogonal to the real quaternion.
    """
    real_part = np.copy(real_part)
    dual_part = np.copy(dual_part)

    real_norm = np.linalg.norm(real_part, axis=1)

    # special case: real quaternion is 0, we set it to identity
    zero_real_mask = real_norm == 0.0
    real_part[zero_real_mask, :4] = [0., 0., 0., 1.]
    real_norm[zero_real_mask] = 1.0

    # 1. ensure unit real quaternion
    real_part /= real_norm[:, np.newaxis]
    dual_part /= real_norm[:, np.newaxis]

    # 2. ensure orthogonality of real and dual quaternion
    dual_part -= np.sum(real_part * dual_part, axis=1)[:, np.newaxis] * real_part

    return real_part, dual_part


def _create_transformation_matrix(translations, rotation_matrices, single):
    """Create a matrix from translations and rotation matrices.

    Parameters
    ----------
    translations : array_like, shape (N, 3) or (3,)
        A stack of translation vectors.
    rotation_matrices : array_like, shape (N, 3, 3) or (3, 3)
        A stack of rotation matrices.
    single : bool
        Whether the output should be a single matrix or a stack of matrices.

    Returns
    -------
    matrix : numpy.ndarray, shape (N, 4, 4)
        A stack of transformation matrices.
    """
    if translations.ndim == 1:
        translations = translations[np.newaxis, ...]
    if rotation_matrices.ndim == 2:
        rotation_matrices = rotation_matrices[np.newaxis, ...]
    if not len(rotation_matrices) == len(translations):
        raise ValueError("The number of rotation matrices and translations must "
                         "be the same.")

    matrix = np.empty((len(translations), 4, 4), dtype=float)
    matrix[:, :3, :3] = rotation_matrices
    matrix[:, :3, 3] = translations
    matrix[:, 3, :3] = 0.0
    matrix[:, 3, 3] = 1.0
    if single:
        return matrix[0]
    else:
        return matrix


cdef class RigidTransform:
    cdef double[:, :, :] _matrix
    cdef bint _single

    def __init__(self, matrix, normalize=True, copy=True):
        self._single = False
        matrix = np.asarray(matrix, dtype=float)

        if (matrix.ndim not in [2, 3]
                or matrix.shape[-1] != 4
                or matrix.shape[-2] != 4):
            raise ValueError("Expected `matrix` to have shape (4, 4), or (N, 4, 4), "
                             f"got {matrix.shape}.")

        # If a single transformation matrix is given, convert it to a
        # 2D 1 x 4 x 4 matrix but set self._single to True so that we can
        # return appropriate objects in the `to_...` methods
        if matrix.shape == (4, 4):
            matrix = matrix[None, :, :]
            self._single = True

        if normalize or copy:
            matrix = matrix.copy()

        # Rigid transforms have the following matrix representation:
        # [R | t]
        # [0 | 1]
        # where R is a 3x3 orthonormal rotation matrix and t is a 3x1 translation
        # vector. The last row is always [0, 0, 0, 1] exactly

        # Check the last row. Exact match is required, as this last row should
        # not accumulate floating point error during composition.
        last_row_ok = np.all(matrix[:, 3, :] == np.array([0, 0, 0, 1]), axis=1)
        if np.any(~last_row_ok):
            ind = np.where(~last_row_ok)[0][0]
            raise ValueError(f"Expected last row of transformation matrix {ind} to be "
                             f"exactly [0, 0, 0, 1], got {matrix[ind, 3, :]}.")

        # The Rotation.from_matrix() method orthogonalizes the rotation
        # component of the transformation matrix. While this does have some
        # overhead in converting a rotation matrix to a quaternion and back, it
        # allows for skipping singular value decomposition for near-orthogonal
        # matrices, which is a computationally expensive operation.
        if normalize:
            matrix[:, :3, :3] = Rotation.from_matrix(matrix[:, :3, :3]).as_matrix()

        self._matrix = matrix

    def __repr__(self):
        m = f"{self.as_matrix()!r}".splitlines()
        # bump indent (+27 characters)
        m[1:] = [" " * 27 + m[i] for i in range(1, len(m))]
        return "RigidTransform.from_matrix(" + "\n".join(m) + ")"

    @cython.embedsignature(True)
    @classmethod
    def from_matrix(cls, matrix):
        return cls(matrix, normalize=True, copy=True)

    @cython.embedsignature(True)
    @classmethod
    def from_rotation(cls, rotation):
        if not isinstance(rotation, Rotation):
            raise ValueError("Expected `rotation` to be a `Rotation` instance, "
                             f"got {type(rotation)}.")
        rotmat = rotation.as_matrix()
        if rotation.single:
            rotmat = rotmat[None, :, :]
        num_transforms = len(rotmat)
        matrix = np.zeros((num_transforms, 4, 4), dtype=float)
        matrix[:, :3, :3] = rotmat
        matrix[:, 3, 3] = 1
        if rotation.single:
            matrix = matrix[0]
        return cls(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_translation(cls, translation):
        translation = np.asarray(translation, dtype=float)

        if translation.ndim not in [1, 2] or translation.shape[-1] != 3:
            raise ValueError("Expected `translation` to have shape (3,), or (N, 3), "
                             f"got {translation.shape}.")

        # If a single translation vector is given, convert it to a 2D 1 x 3 matrix
        single = False
        if translation.shape == (3,):
            translation = translation[None, :]
            single = True

        num_translations = translation.shape[0]

        matrix = np.repeat(np.eye(4, dtype=float)[None, :, :], num_translations, axis=0)
        matrix[:, :3, 3] = translation

        if single:
            matrix = matrix[0]
        return cls(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_components(cls, translation, rotation):
        return cls.from_translation(translation) * cls.from_rotation(rotation)

    @cython.embedsignature(True)
    @classmethod
    def from_exp_coords(cls, exp_coords):
        exp_coords = np.asarray(exp_coords, dtype=float)

        if exp_coords.ndim not in [1, 2] or exp_coords.shape[-1] != 6:
            raise ValueError(
                "Expected `exp_coords` to have shape (6,), or (N, 6), "
                f"got {exp_coords.shape}.")

        single = exp_coords.ndim == 1
        exp_coords = np.atleast_2d(exp_coords)

        rot_vec = exp_coords[:, :3]
        rotations = Rotation.from_rotvec(rot_vec)

        translations = np.einsum('ijk,ik->ij',
                                 _compute_se3_exp_translation_transform(rot_vec),
                                 exp_coords[:, 3:])
        matrix = _create_transformation_matrix(translations, rotations.as_matrix(),
                                               single)

        return cls(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def from_dual_quat(cls, dual_quat, *, scalar_first=False):
        dual_quat = np.asarray(dual_quat, dtype=float)

        if dual_quat.ndim not in [1, 2] or dual_quat.shape[-1] != 8:
            raise ValueError(
                "Expected `dual_quat` to have shape (8,), or (N, 8), "
                f"got {dual_quat.shape}.")

        single = dual_quat.ndim == 1
        dual_quat = np.atleast_2d(dual_quat)

        real_part = dual_quat[:, :4]
        dual_part = dual_quat[:, 4:]
        if scalar_first:
            real_part = np.roll(real_part, -1, axis=1)
            dual_part = np.roll(dual_part, -1, axis=1)

        real_part, dual_part = _normalize_dual_quaternion(real_part, dual_part)

        rotation = Rotation.from_quat(real_part)
        translation = 2.0 * np.asarray(
            compose_quat(dual_part, rotation.inv().as_quat()))[:, :3]
        matrix = _create_transformation_matrix(translation, rotation.as_matrix(),
                                               single)

        return cls(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    @classmethod
    def identity(cls, num=None):
        if num is None:  # single
            return cls(np.eye(4, dtype=float), normalize=False, copy=False)
        else:
            return cls(np.eye(4, dtype=float)[None, :, :].repeat(num, axis=0),
                       normalize=False, copy=False)


    @cython.embedsignature(True)
    @classmethod
    def concatenate(cls, transforms):
        if isinstance(transforms, RigidTransform):
            return cls(transforms._matrix, normalize=False, copy=True)

        if not all(isinstance(x, RigidTransform) for x in transforms):
            raise TypeError("input must contain RigidTransform objects only")

        ms = [x.as_matrix()[np.newaxis, :, :] if x.single else x.as_matrix()
              for x in transforms]
        return cls(np.concatenate(ms), normalize=False, copy=False)

    @cython.embedsignature(True)
    def as_matrix(self):
        if self._single:
            return np.array(self._matrix[0])
        else:
            return np.array(self._matrix)

    @cython.embedsignature(True)
    def as_components(self):
        return self.translation, self.rotation

    @cython.embedsignature(True)
    def as_exp_coords(self):
        exp_coords = np.empty((self._matrix.shape[0], 6), dtype=float)
        rotations = Rotation.from_matrix(self._matrix[:, :3, :3])
        rot_vec = rotations.as_rotvec()
        exp_coords[:, :3] = rot_vec
        exp_coords[:, 3:] = np.einsum('ijk,ik->ij',
                                      _compute_se3_log_translation_transform(rot_vec),
                                      self._matrix[:, :3, 3])

        if self._single:
            return exp_coords[0]
        else:
            return exp_coords

    @cython.embedsignature(True)
    def as_dual_quat(self, *, scalar_first=False):
        rotations = Rotation.from_matrix(self._matrix[:, :3, :3])
        real_parts = rotations.as_quat()

        pure_translation_quats = np.empty((len(self._matrix), 4), dtype=float)
        pure_translation_quats[:, :3] = self._matrix[:, :3, 3]
        pure_translation_quats[:, 3] = 0.0

        dual_parts = 0.5 * np.asarray(
            compose_quat(pure_translation_quats, real_parts))

        if scalar_first:
            real_parts = np.roll(real_parts, 1, axis=1)
            dual_parts = np.roll(dual_parts, 1, axis=1)

        dual_quats = np.hstack((real_parts, dual_parts))

        if self._single:
            return dual_quats[0]
        else:
            return dual_quats

    @cython.embedsignature(True)
    def __len__(self):
        if self._single:
            raise TypeError("Single transform has no len().")

        return self._matrix.shape[0]

    @cython.embedsignature(True)
    def __getitem__(self, indexer):
        if self._single:
            raise TypeError("Single transform is not subscriptable.")

        # Convert memoryview to numpy array before indexing:
        arr = np.asarray(self._matrix)
        return self.__class__(arr[indexer], copy=False)

    @cython.embedsignature(True)
    def __setitem__(self, indexer, value):
        if self._single:
            raise TypeError("Single transform is not subscriptable.")

        if not isinstance(value, self.__class__):
            raise TypeError("value must be a RigidTransform object")

        arr = np.asarray(self._matrix)
        arr[indexer] = value.as_matrix()
        self._matrix = arr

    @cython.embedsignature(True)
    def __mul__(RigidTransform self, RigidTransform other):
        len_self = len(self._matrix)
        len_other = len(other._matrix)
        if not(len_self == 1 or len_other == 1 or len_self == len_other):
            raise ValueError("Expected equal number of transforms in both or a "
                             "single transform in either object, got "
                             f"{len_self} transforms in first and {len_other}"
                             "transforms in second object.")

        result = np.matmul(self._matrix, other._matrix)
        if self._single and other._single:
            result = result[0]
        return self.__class__(result, copy=False)

    @cython.embedsignature(True)
    def __pow__(RigidTransform self, float n):
        # Exact short-cuts
        if n == 0:
            return self.__class__.identity(None if self._single else len(self._matrix))
        elif n == -1:
            return self.inv()
        elif n == 1:
            if self._single:
                return self.__class__(self._matrix[0], copy=True)
            else:
                return self.__class__(self._matrix, copy=True)
        else:
            return self.__class__.from_exp_coords(self.as_exp_coords() * n)

    @cython.embedsignature(True)
    def inv(self):
        r_inv = np.swapaxes(self._matrix[:, :3, :3], 1, 2)  # Transpose to invert
        # This einsum performs element-wise matrix multiplication
        t_inv = -np.einsum('ijk,ik->ij', r_inv, self._matrix[:, :3, 3])

        matrix = _create_transformation_matrix(t_inv, r_inv, self._single)
        return self.__class__(matrix, normalize=False, copy=False)

    @cython.embedsignature(True)
    def apply(self, vector, inverse=False):
        vector = np.asarray(vector, dtype=float)
        if vector.ndim not in [1, 2] or vector.shape[-1] != 3:
            raise ValueError("Expected vector to have shape (N, 3), or (3,), "
                             f"got {vector.shape}.")

        vector_single = False
        if vector.ndim == 1:
            vector = vector[None, :]
            vector_single = True

        vector = np.hstack([vector, np.ones((vector.shape[0], 1))])

        if inverse:
            m = self.inv().as_matrix()
            if self._single:
                m = m[None, :, :]
        else:
            m = self._matrix

        # This einsum performs matrix multiplication of each of the (4, 4)
        # matrices in `m` with the (4,) vectors in `vector`, with proper
        # broadcasting for different dimensions of `m` and `vector`.
        res = np.einsum('ijk,ik->ij', m, vector)[:, :3]
        if self._single and vector_single:
            return res[0]
        else:
            return res

    @cython.embedsignature(True)
    @property
    def rotation(self):
        if self._single:
            return Rotation.from_matrix(self._matrix[0, :3, :3])
        else:
            return Rotation.from_matrix(self._matrix[:, :3, :3])

    @cython.embedsignature(True)
    @property
    def translation(self):
        if self._single:
            return np.array(self._matrix[0, :3, 3])
        else:
            return np.array(self._matrix[:, :3, 3])

    @cython.embedsignature(True)
    @property
    def single(self):
        return self._single
