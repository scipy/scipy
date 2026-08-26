import sys
import pytest
import warnings
# importing `scipy.stats` does not emit warning
with warnings.catch_warnings():
    warnings.simplefilter("error")
    from scipy import stats
    import scipy.stats  # noqa: F401


class TestMstatsDeprecation:
    def forget_mstats(self):
        sys.modules.pop("scipy.stats.mstats", None)
        stats.__dict__.pop("mstats", None)
        globals().pop("mstats", None)

    def test_import_message(self):
        message = "`scipy.stats.mstats` is deprecated..."

        self.forget_mstats()
        with pytest.warns(DeprecationWarning, match=message):
            stats.mstats

        self.forget_mstats()
        with pytest.warns(DeprecationWarning, match=message):
            import scipy.stats.mstats  # noqa: F401

        self.forget_mstats()
        with pytest.warns(DeprecationWarning, match=message):
            from scipy.stats import mstats  # noqa: F401

    def test_function_warnings(self):
        # Spot check of the various branches of _get_deprecation_message
        x = [1., 2., 3.]

        message = "SciPy offers no replacement for this function."
        with pytest.warns(DeprecationWarning, match=message):
            stats.mstats.argstoarray([x])

        message = "use `scipy.stats.gmean` with MArray\\(s\\)..."
        with pytest.warns(DeprecationWarning, match=message):
            stats.mstats.gmean(x)

        message = "use `scipy.stats.quantile` with MArray\\(s\\)..."
        with pytest.warns(DeprecationWarning, match=message):
            stats.mstats.hdquantiles(x)

        message = "use `scipy.stats.obrientransform` with regular NumPy array\\(s\\)..."
        with pytest.warns(DeprecationWarning, match=message):
            stats.mstats.obrientransform([x])

        message = "See the Trimming and winsorization transition guide..."
        with pytest.warns(DeprecationWarning, match=message):
            stats.mstats.trima([x])
