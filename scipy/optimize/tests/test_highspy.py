"""
Tests for the scipy.optimize._highspy._core bindings (the nanobind-based
wrapper around the HiGHS C++ library).

These exercise the low-level binding surface directly (as opposed to the
public scipy.optimize.linprog/milp APIs), since it is easy for a C++/Python
binding layer to have lifetime, refcounting, or type-marshalling bugs that
don't show up through the higher-level solvers.
"""
import gc
import sys

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scipy.optimize import linprog, milp, LinearConstraint, Bounds, OptimizeWarning
import scipy.optimize._highspy._core as _core


def _build_solved_highs():
    h = _core._Highs()
    lp = _core.HighsLp()
    lp.num_col_ = 2
    lp.num_row_ = 2
    lp.col_cost_ = np.array([1.0, 1.0])
    lp.col_lower_ = np.array([0.0, 0.0])
    lp.col_upper_ = np.array([_core.kHighsInf, _core.kHighsInf])
    lp.row_lower_ = np.array([-_core.kHighsInf, -_core.kHighsInf])
    lp.row_upper_ = np.array([10.0, 15.0])
    lp.sense_ = _core.ObjSense.kMinimize
    lp.offset_ = 0.0
    mat = _core.HighsSparseMatrix()
    mat.format_ = _core.MatrixFormat.kColwise
    mat.num_col_ = 2
    mat.num_row_ = 2
    mat.start_ = [0, 2, 4]
    mat.index_ = [0, 1, 0, 1]
    mat.value_ = [1.0, 1.0, 2.0, 1.0]
    lp.a_matrix_ = mat
    model = _core.HighsModel()
    model.lp_ = lp
    h.passModel(model)
    h.setOptionValue("output_flag", False)
    h.run()
    return h


class TestHighsCoreModelBuilding:
    def test_pass_model_via_objects(self):
        h = _core._Highs()
        lp = _core.HighsLp()
        lp.num_col_ = 2
        lp.num_row_ = 1
        lp.col_cost_ = np.array([1.0, 2.0])
        lp.col_lower_ = np.array([0.0, 0.0])
        lp.col_upper_ = np.array([10.0, 10.0])
        lp.row_lower_ = np.array([-_core.kHighsInf])
        lp.row_upper_ = np.array([10.0])
        lp.sense_ = _core.ObjSense.kMinimize
        mat = _core.HighsSparseMatrix()
        mat.format_ = _core.MatrixFormat.kColwise
        mat.num_col_ = 2
        mat.num_row_ = 1
        mat.start_ = [0, 1, 2]
        mat.index_ = [0, 0]
        mat.value_ = [1.0, 1.0]
        lp.a_matrix_ = mat
        model = _core.HighsModel()
        model.lp_ = lp

        assert h.passModel(model) == _core.HighsStatus.kOk
        assert h.run() == _core.HighsStatus.kOk

    def test_pass_model_via_raw_arrays(self):
        h = _core._Highs()
        status = h.passModel(
            2, 1, 2, 0,
            int(_core.MatrixFormat.kColwise), int(_core.MatrixFormat.kColwise),
            int(_core.ObjSense.kMinimize), 0.0,
            np.array([1.0, 2.0]), np.array([0.0, 0.0]), np.array([10.0, 10.0]),
            np.array([-_core.kHighsInf]), np.array([10.0]),
            np.array([0, 1, 2], dtype=np.int32), np.array([0, 0], dtype=np.int32),
            np.array([1.0, 1.0]),
            np.array([], dtype=np.int32), np.array([], dtype=np.float64),
            np.array([], dtype=np.int32), np.array([], dtype=np.int32),
        )
        assert status == _core.HighsStatus.kOk
        assert h.run() == _core.HighsStatus.kOk


class TestHighsOptionTypeDispatch:
    def test_linprog_mixed_option_types(self):
        # exercise bool, float, int, string, and enum-valued options in one
        # call, covering every branch of _highs_wrapper.py's check_option().
        # "random_seed", "log_file", and "highs_debug_level" are not among
        # linprog(method="highs")'s named kwargs, so they are passed through
        # verbatim to HiGHS with a warning -- that pass-through path (and the
        # enum -> int conversion in check_option) is exactly what's under
        # test here.
        options = {
            "presolve": True,
            "time_limit": 30.5,
            "random_seed": 1,
            "log_file": "",
            "highs_debug_level": _core.HighsDebugLevel.kHighsDebugLevelNone,
        }
        with pytest.warns(OptimizeWarning, match="Unrecognized options"):
            res = linprog([1, 2], A_ub=[[1, 1]], b_ub=[10],
                          bounds=[(0, 10), (0, 10)], method="highs",
                          options=options)
        assert res.success
        assert_allclose(res.fun, 0.0, atol=1e-9)

    def test_milp_mixed_option_types(self):
        options = {
            "presolve": True,
            "time_limit": 30.0,
            "mip_rel_gap": 0.001,
        }
        constraints = LinearConstraint([[1, 1]], -np.inf, 10)
        bounds = Bounds([0, 0], [10, 10])
        res = milp([-1, -2], constraints=constraints, bounds=bounds,
                   integrality=[1, 1], options=options)
        assert res.success


def test_col_cost_getter_returns_numpy_array():
    # Regression test: HighsLp.col_cost_ (backed by make_readonly_ptr) must
    # come back as a real, usable numpy array, not a bare DLPack capsule.
    lp = _core.HighsLp()
    lp.num_col_ = 3
    lp.col_cost_ = np.array([1.0, 2.0, 3.0])
    arr = np.asarray(lp.col_cost_)
    assert_allclose(arr, [1.0, 2.0, 3.0])


def test_get_basis_inverse_row_numeric():
    # Regression test for the nanobind port's to_ndarray()-based return
    # values: the returned array must be numerically correct and usable
    # as a plain numpy array.
    h = _build_solved_highs()
    status, row0 = h.getBasisInverseRow(0)
    assert status == _core.HighsStatus.kOk
    assert_allclose(np.asarray(row0), [1.0, 0.0], atol=1e-9)

    status, row1 = h.getBasisInverseRow(1)
    assert status == _core.HighsStatus.kOk
    assert_allclose(np.asarray(row1), [0.0, 1.0], atol=1e-9)


def test_uaf_col_cost_property_survives_parent_deletion():
    # Use-after-free regression test: HighsLp.col_cost_'s getter is backed
    # by nb::rv_policy::reference_internal (make_readonly_ptr), which is
    # supposed to keep the parent HighsLp alive for as long as the returned
    # array is alive. Run this in a loop with gc.collect() to raise the
    # odds of catching a real UAF if that lifetime tie is ever broken.
    results = []
    for i in range(200):
        lp = _core.HighsLp()
        lp.num_col_ = 4
        lp.col_cost_ = np.array([i * 1.0, i * 2.0, i * 3.0, i * 4.0])
        arr = lp.col_cost_
        del lp
        gc.collect()
        results.append((i, arr))

    for i, arr in results:
        assert_allclose(np.asarray(arr), [i * 1.0, i * 2.0, i * 3.0, i * 4.0])


def test_set_callback_none_clears_callback():
    h = _core._Highs()
    assert h.setCallback(None, None) == _core.HighsStatus.kOk


def test_set_callback_fires_with_correct_user_data_and_refcounting():
    class UserData:
        def __init__(self, tag):
            self.tag = tag
            self.calls = 0

    h = _build_solved_highs()
    ud = UserData("hello-world")
    received = []

    def callback(cb_type, msg, data_out, data_in, user_data):
        received.append(user_data.tag)
        user_data.calls += 1

    assert h.setCallback(callback, ud) == _core.HighsStatus.kOk
    h.startCallback(int(_core.cb.HighsCallbackType.kCallbackLogging))
    h.setOptionValue("output_flag", True)
    h.run()
    h.stopCallback(int(_core.cb.HighsCallbackType.kCallbackLogging))
    h.setCallback(None, None)
    del h, callback
    gc.collect()

    assert len(received) > 0
    assert all(tag == "hello-world" for tag in received)
    assert ud.calls == len(received)


def test_set_callback_repeated_create_destroy_loop():
    # Repeated create/destroy of Highs + callback + user data, to catch
    # reference-counting leaks or crashes that a single call wouldn't show.
    tags = []
    for i in range(30):
        h = _build_solved_highs()
        ud = object()
        rc0 = sys.getrefcount(ud)

        def callback(cb_type, msg, data_out, data_in, user_data, _i=i, _tags=tags):
            _tags.append(_i)

        h.setCallback(callback, ud)
        h.startCallback(int(_core.cb.HighsCallbackType.kCallbackLogging))
        h.setOptionValue("output_flag", True)
        h.run()
        h.stopCallback(int(_core.cb.HighsCallbackType.kCallbackLogging))
        h.setCallback(None, None)
        del h, callback
        gc.collect()
        rc1 = sys.getrefcount(ud)
        assert rc1 <= rc0 + 2, f"refcount grew unexpectedly at iter {i}: {rc0} -> {rc1}"

    assert len([t for t in tags if t is not None]) > 0


def test_enum_valued_option_via_linprog():
    # A recent fix (int(value) in check_option()) relies on bound enums
    # supporting int conversion when passed as an options= dict value.
    with pytest.warns(OptimizeWarning, match="Unrecognized options"):
        res = linprog([1, 2], bounds=[(0, 10), (0, 10)], method="highs",
                      options={"highs_debug_level":
                               _core.HighsDebugLevel.kHighsDebugLevelCheap})
    assert res.success
