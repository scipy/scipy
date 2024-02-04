By Daniel Lemire

See https://github.com/fastfloat/fast_float


# fast_float number parsing library: 4x faster than strtod

The fast_float library provides fast header-only implementations for the C++ from_chars
functions for `float` and `double` types.  These functions convert ASCII strings representing
decimal values (e.g., `1.3e10`) into binary types. We provide exact rounding (including
round to even). In our experience, these `fast_float` functions many times faster than comparable number-parsing functions from existing C++ standard libraries.

Specifically, `fast_float` provides the following two functions with a C++17-like syntax (the library itself only requires C++11):

```C++
from_chars_result from_chars(const char* first, const char* last, float& value, ...);
from_chars_result from_chars(const char* first, const char* last, double& value, ...);
```

The return type (`from_chars_result`) is defined as the struct:
```C++
struct from_chars_result {
    const char* ptr;
    std::errc ec;
};
```

It parses the character sequence [first,last) for a number. It parses floating-point numbers expecting
a locale-independent format equivalent to the C++17 from_chars function.
The resulting floating-point value is the closest floating-point values (using either float or double),
using the "round to even" convention for values that would otherwise fall right in-between two values.
That is, we provide exact parsing according to the IEEE standard.


Given a successful parse, the pointer (`ptr`) in the returned value is set to point right after the
parsed number, and the `value` referenced is set to the parsed value. In case of error, the returned
`ec` contains a representative error, otherwise the default (`std::errc()`) value is stored.

The implementation does not throw and does not allocate memory (e.g., with `new` or `malloc`).

It will parse infinity and nan values.