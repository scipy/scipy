# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 10:54:01 2021

@author: matth
"""

from scipy import stats
from scipy.stats.tests import common_tests
# print(stats.norm.moment(n=1))  # 0.0; want DeprecationWarning
# print(stats.t.moment(n=1, df=10, loc=1))  # 1.0; want DeprecationWarning


common_tests.check_deprecation_warning_gh5982(stats.norm, tuple(), 'norm')

# using n as a keyword argument for moment has NOT been possible, so there
# should be no deprecation warnings
# print(stats.binom.moment(1, 10, 0.5))  # 5.0
# print(stats.binom.moment(1, n=10, p=0.5))  # TypeError: moment() got multiple values for argument 'n'
# print(stats.binom.moment(n=10, p=0.5))  # TypeError: moment() missing...

# print(stats.binom.moment(n=1, n=10, p=0.5))  # SyntaxError: keyword argument repeated: n
# print(stats.binom.moment(order=1, n=10, p=0.5))  # this _should_ work

# print(stats.norm.interval(alpha=0.5))  # (-0.6744897501960817, 0.6744897501960817)
# print(stats.t.interval(alpha=0.5, df=10, loc=1))  # (0.30018793868757077, 1.6998120613124292)

# a, b = 1.8, -0.5
# print(stats.levy_stable.interval(0.5, a, b))  # (-0.8852766959515754, 1.0409855813082687)
# print(stats.levy_stable.interval(0.5, alpha=a, beta=b))  # previously TypeError: interval() got multiple values for argument 'alpha'
# # print(stats.levy_stable.interval(alpha=a, beta=b))  # TypeError: interval() missing 1 required positional argument: `confidence`
# print(stats.levy_stable.interval(alpha=0.5, alpha=a, beta=b))  # SyntaxError: keyword argument repeated: alpha
# print(stats.levy_stable.interval(confidence=0.5, alpha=a, beta=b))  # SyntaxError: keyword argument repeated: alpha


# moment
# n as moment kwd, n as shape kwd  # SyntaxError (keyword argument repeated) in master and here
# n as moment kwd, n as shape pos  # SyntaxError (positional argument follows keyword argument) in master and here
# order as moment kwd, n as shape kwd  # works fine
# order as moment kwd, n as shape pos  # SyntaxError (positional argument follows keyword argument) in master and here
# first argument moment pos, n as shape kwd  # previously TypeError: moment() got multiple values for argument 'n', now works fine
# first argument moment pos, n as shape pos  # always worked fine