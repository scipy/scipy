.. _writing-test-tips:

===========================
Tips for writing tests
===========================

Assert function selection guideline
------------------------------------
When writing unit tests, the guideline for selecting which assert function
to use are as follows:

- ``xp_assert_close``: for floating-point comparison

- ``xp_assert_equal``: for integer arrays

- ``assert`` for anything else (boolean check, etc.)