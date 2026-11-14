#pragma once

#include <cassert>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>

/**
 * This namespace defines wrappers for std::ifstream, std::ofstream, and
 * std::fstream objects. The wrappers perform the following steps:
 * - check the open modes make sense
 * - check that the call to open() is successful
 * - (for input streams) check that the opened file is peek-able
 * - turn on the badbit in the exception mask
 */
namespace strict_fstream
{

// Help people out a bit, it seems like this is a common recommendation since
// musl breaks all over the place.
#if defined(__NEED_size_t) && !defined(__MUSL__)
#warning "It seems to be recommended to patch in a define for __MUSL__ if you use musl globally: https://www.openwall.com/lists/musl/2013/02/10/5"
#define __MUSL__
#endif

// Workaround for broken musl implementation
// Since musl insists that they are perfectly compatible, ironically enough,
// they don't officially have a __musl__ or similar. But __NEED_size_t is defined in their
// relevant header (and not in working implementations), so we can use that.
#ifdef __MUSL__
#warning "Working around broken strerror_r() implementation in musl, remove when musl is fixed"
#endif

// Non-gnu variants of strerror_* don't necessarily null-terminate if
// truncating, so we have to do things manually.
inline std::string trim_to_null(const std::vector<char> &buff)
{
    std::string ret(buff.begin(), buff.end());

    const std::string::size_type pos = ret.find('\0');
    if (pos == std::string::npos) {
        ret += " [...]"; // it has been truncated
    } else {
        ret.resize(pos);
    }
    return ret;
}

/// Overload of error-reporting function, to enable use with VS and non-GNU
/// POSIX libc's
/// Ref:
///   - http://stackoverflow.com/a/901316/717706
static std::string strerror()
{
    // Can't use std::string since we're pre-C++17
    std::vector<char> buff(256, '\0');

#ifdef _WIN32
    // Since strerror_s might set errno itself, we need to store it.
    const int err_num = errno;
    if (strerror_s(buff.data(), buff.size(), err_num) != 0) {
        return trim_to_null(buff);
    } else {
        return "Unknown error (" + std::to_string(err_num) + ")";
    }
#elif ((_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600 || defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__)) && ! _GNU_SOURCE) || defined(__MUSL__)
// XSI-compliant strerror_r()
    const int err_num = errno; // See above
    if (strerror_r(err_num, buff.data(), buff.size()) == 0) {
        return trim_to_null(buff);
    } else {
        return "Unknown error (" + std::to_string(err_num) + ")";
    }
#else
// GNU-specific strerror_r()
    char * p = strerror_r(errno, &buff[0], buff.size());
    return std::string(p, std::strlen(p));
#endif
}

/// Exception class thrown by failed operations.
class Exception
    : public std::exception
{
public:
    Exception(const std::string& msg) : _msg(msg) {}
    const char * what() const noexcept { return _msg.c_str(); }
private:
    std::string _msg;
}; // class Exception

namespace detail
{

struct static_method_holder
{
    static std::string mode_to_string(std::ios_base::openmode mode)
    {
        static const int n_modes = 6;
        static const std::ios_base::openmode mode_val_v[n_modes] =
            {
                std::ios_base::in,
                std::ios_base::out,
                std::ios_base::app,
                std::ios_base::ate,
                std::ios_base::trunc,
                std::ios_base::binary
            };

        static const char * mode_name_v[n_modes] =
            {
                "in",
                "out",
                "app",
                "ate",
                "trunc",
                "binary"
            };
        std::string res;
        for (int i = 0; i < n_modes; ++i)
        {
            if (mode & mode_val_v[i])
            {
                res += (! res.empty()? "|" : "");
                res += mode_name_v[i];
            }
        }
        if (res.empty()) res = "none";
        return res;
    }
    static void check_mode(const std::string& filename, std::ios_base::openmode mode)
    {
        if ((mode & std::ios_base::trunc) && ! (mode & std::ios_base::out))
        {
            throw Exception(std::string("strict_fstream: open('") + filename + "'): mode error: trunc and not out");
        }
        else if ((mode & std::ios_base::app) && ! (mode & std::ios_base::out))
        {
            throw Exception(std::string("strict_fstream: open('") + filename + "'): mode error: app and not out");
        }
        else if ((mode & std::ios_base::trunc) && (mode & std::ios_base::app))
        {
            throw Exception(std::string("strict_fstream: open('") + filename + "'): mode error: trunc and app");
        }
     }
    static void check_open(std::ios * s_p, const std::string& filename, std::ios_base::openmode mode)
    {
        if (s_p->fail())
        {
            throw Exception(std::string("strict_fstream: open('")
                            + filename + "'," + mode_to_string(mode) + "): open failed: "
                            + strerror());
        }
    }
    static void check_peek(std::istream * is_p, const std::string& filename, std::ios_base::openmode mode)
    {
        bool peek_failed = true;
        try
        {
            is_p->peek();
            peek_failed = is_p->fail();
        }
        catch (const std::ios_base::failure &) {}
        if (peek_failed)
        {
            throw Exception(std::string("strict_fstream: open('")
                            + filename + "'," + mode_to_string(mode) + "): peek failed: "
                            + strerror());
        }
        is_p->clear();
    }
}; // struct static_method_holder

} // namespace detail

class ifstream
    : public std::ifstream
{
public:
    ifstream() = default;
    ifstream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in)
    {
        open(filename, mode);
    }
    void open(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in)
    {
        mode |= std::ios_base::in;
        exceptions(std::ios_base::badbit);
        detail::static_method_holder::check_mode(filename, mode);
        std::ifstream::open(filename, mode);
        detail::static_method_holder::check_open(this, filename, mode);
        detail::static_method_holder::check_peek(this, filename, mode);
    }
}; // class ifstream

class ofstream
    : public std::ofstream
{
public:
    ofstream() = default;
    ofstream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out)
    {
        open(filename, mode);
    }
    void open(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out)
    {
        mode |= std::ios_base::out;
        exceptions(std::ios_base::badbit);
        detail::static_method_holder::check_mode(filename, mode);
        std::ofstream::open(filename, mode);
        detail::static_method_holder::check_open(this, filename, mode);
    }
}; // class ofstream

class fstream
    : public std::fstream
{
public:
    fstream() = default;
    fstream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in)
    {
        open(filename, mode);
    }
    void open(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in)
    {
        if (! (mode & std::ios_base::out)) mode |= std::ios_base::in;
        exceptions(std::ios_base::badbit);
        detail::static_method_holder::check_mode(filename, mode);
        std::fstream::open(filename, mode);
        detail::static_method_holder::check_open(this, filename, mode);
        detail::static_method_holder::check_peek(this, filename, mode);
    }
}; // class fstream

} // namespace strict_fstream

