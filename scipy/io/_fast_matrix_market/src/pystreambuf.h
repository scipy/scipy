// SPDX-License-Identifier: BSD-3-Clause
/*
Based on https://gist.github.com/asford/544323a5da7dddad2c9174490eb5ed06

Original license text
---------------------

This component utilizes components derived from cctbx, available at
http://cci.lbl.gov/cctbx_sources/boost_adaptbx/python_streambuf.h

*** License agreement ***

cctbx Copyright (c) 2006, The Regents of the University of
California, through Lawrence Berkeley National Laboratory (subject to
receipt of any required approvals from the U.S. Dept. of Energy).  All
rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes,
patches, or upgrades to the features, functionality or performance of
the source code ("Enhancements") to anyone; however, if you choose to
make your Enhancements available either publicly, or directly to
Lawrence Berkeley National Laboratory, without imposing a separate
written license agreement for such Enhancements, then you hereby grant
the following license: a  non-exclusive, royalty-free perpetual license
to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.

*/

#pragma once

#include <pybind11/pybind11.h>

#include <streambuf>
#include <iostream>

namespace py = pybind11;

/// A stream buffer getting data from and putting data into a Python file object
/** The aims are as follow:

    - Given a C++ function acting on a standard stream, e.g.

      \code
      void read_inputs(std::istream& input) {
        ...
        input >> something >> something_else;
      }
      \endcode

      and given a piece of Python code which creates a file-like object,
      to be able to pass this file object to that C++ function, e.g.

      \code
      import gzip
      gzip_file_obj = gzip.GzipFile(...)
      read_inputs(gzip_file_obj)
      \endcode

      and have the standard stream pull data from and put data into the Python
      file object.

    - When Python \c read_inputs() returns, the Python object is able to
      continue reading or writing where the C++ code left off.

    - Operations in C++ on mere files should be competitively fast compared
      to the direct use of \c std::fstream.


    \b Motivation

      - the standard Python library offer of file-like objects (files,
        compressed files and archives, network, ...) is far superior to the
        offer of streams in the C++ standard library and Boost C++ libraries.

      - i/o code involves a fair amount of text processing which is more
        efficiently prototyped in Python but then one may need to rewrite
        a time-critical part in C++, in as seamless a manner as possible.

    \b Usage

    This is 2-step:

      - a trivial wrapper function

        \code
          using boost_adaptbx::python::streambuf;
          void read_inputs_wrapper(streambuf& input)
          {
            streambuf::istream is(input);
            read_inputs(is);
          }

          def("read_inputs", read_inputs_wrapper);
        \endcode

        which has to be written every time one wants a Python binding for
        such a C++ function.

      - the Python side

        \code
          from boost.python import streambuf
          read_inputs(streambuf(python_file_obj=obj, buffer_size=1024))
        \endcode

        \c buffer_size is optional. See also: \c default_buffer_size

  Note: references are to the C++ standard (the numbers between parentheses
  at the end of references are margin markers).
*/

namespace pystream{

class streambuf : public std::basic_streambuf<char>
{
  private:
    typedef std::basic_streambuf<char> base_t;

  public:
    /* The syntax
        using base_t::char_type;
       would be nicer but Visual Studio C++ 8 chokes on it
    */
    typedef base_t::char_type   char_type;
    typedef base_t::int_type    int_type;
    typedef base_t::pos_type    pos_type;
    typedef base_t::off_type    off_type;
    typedef base_t::traits_type traits_type;

    /// The default size of the read and write buffer.
    /** They are respectively used to buffer data read from and data written to
        the Python file object. It can be modified from Python.
    */
    static inline std::size_t default_buffer_size = 1024;

    /// Construct from a Python file object
    /** if buffer_size is 0 the current default_buffer_size is used.
    */
    streambuf(
      py::object& python_file_obj,
      std::size_t buffer_size_=0)
    :
      py_read (getattr(python_file_obj, "read", py::none())),
      py_write (getattr(python_file_obj, "write", py::none())),
      py_seek (getattr(python_file_obj, "seek", py::none())),
      py_tell (getattr(python_file_obj, "tell", py::none())),
      buffer_size(buffer_size_ != 0 ? buffer_size_ : default_buffer_size),
      write_buffer(0),
      pos_of_read_buffer_end_in_py_file(0),
      pos_of_write_buffer_end_in_py_file(buffer_size),
      farthest_pptr(0)
    {
      assert(buffer_size != 0);
      /* Some Python file objects (e.g. sys.stdout and sys.stdin)
         have non-functional seek and tell. If so, assign None to
         py_tell and py_seek.
       */
      if (!py_tell.is_none()) {
        try {
          py_tell();
        }
        catch (py::error_already_set& err) {
          py_tell = py::none();
          py_seek = py::none();
          err.restore();
          PyErr_Clear();
        }
      }

      if (!py_write.is_none()) {
        // C-like string to make debugging easier
        write_buffer = new char[buffer_size + 1];
        write_buffer[buffer_size] = '\0';
        setp(write_buffer, write_buffer + buffer_size);  // 27.5.2.4.5 (5)
        farthest_pptr = pptr();
      }
      else {
        // The first attempt at output will result in a call to overflow
        setp(0, 0);
      }

      if (!py_tell.is_none()){
        off_type py_pos = py_tell().cast<off_type>();
        pos_of_read_buffer_end_in_py_file = py_pos;
        pos_of_write_buffer_end_in_py_file = py_pos;
      }
    }

    /// Mundane destructor freeing the allocated resources
    virtual ~streambuf() {
      if (write_buffer) delete[] write_buffer;
    }

    /// C.f. C++ standard section 27.5.2.4.3
    /** It is essential to override this virtual function for the stream
        member function readsome to work correctly (c.f. 27.6.1.3, alinea 30)
     */
    virtual std::streamsize showmanyc() {
      int_type const failure = traits_type::eof();
      int_type status = underflow();
      if (status == failure) return -1;
      return egptr() - gptr();
    }

    /// C.f. C++ standard section 27.5.2.4.3
    virtual int_type underflow() {
      int_type const failure = traits_type::eof();
      if (py_read.is_none()) {
        throw std::invalid_argument(
          "That Python file object has no 'read' attribute");
      }
      read_buffer = py_read(buffer_size);
      char *read_buffer_data;
      py::ssize_t py_n_read;
      if (PYBIND11_BYTES_AS_STRING_AND_SIZE(read_buffer.ptr(),
            &read_buffer_data, &py_n_read) == -1) {
        setg(0, 0, 0);
        throw std::invalid_argument(
          "The method 'read' of the Python file object "
          "did not return a string.");
      }
      off_type n_read = (off_type)py_n_read;
      pos_of_read_buffer_end_in_py_file += n_read;
      setg(read_buffer_data, read_buffer_data, read_buffer_data + n_read);
      // ^^^27.5.2.3.1 (4)
      if (n_read == 0) return failure;
      return traits_type::to_int_type(read_buffer_data[0]);
    }

    /// C.f. C++ standard section 27.5.2.4.5
    virtual int_type overflow(int_type c=traits_type::eof()) {
      if (py_write.is_none()) {
        throw std::invalid_argument(
          "That Python file object has no 'write' attribute");
      }
      farthest_pptr = (std::max)(farthest_pptr, pptr());
      off_type n_written = (off_type)(farthest_pptr - pbase());

      // Write all outstanding bytes. Write in chunks to support 32-bit systems.
      // The py::bytes constructor insists that sizeof(second argument) <= sizeof(ssize_t)
      // This is an issue on 32-bit systems where off_type may be long long but ssize_t is an int.
      for (off_type offset = 0; offset < n_written; ) {
        off_type chunk_len = std::min((n_written - offset), (off_type)(1 << 25));
        py::bytes chunk(pbase() + offset, (int)chunk_len);
        py_write(chunk);
        offset += chunk_len;
      }

      if (!traits_type::eq_int_type(c, traits_type::eof())) {
        char cs = traits_type::to_char_type(c);
        py_write(py::bytes(&cs, 1));
        n_written++;
      }
      if (n_written) {
        pos_of_write_buffer_end_in_py_file += n_written;
        setp(pbase(), epptr());
        // ^^^ 27.5.2.4.5 (5)
        farthest_pptr = pptr();
      }
      return traits_type::eq_int_type(
        c, traits_type::eof()) ? traits_type::not_eof(c) : c;
    }

    /// Update the python file to reflect the state of this stream buffer
    /** Empty the write buffer into the Python file object and set the seek
        position of the latter accordingly (C++ standard section 27.5.2.4.2).
        If there is no write buffer or it is empty, but there is a non-empty
        read buffer, set the Python file object seek position to the
        seek position in that read buffer.
    */
    virtual int sync() {
      int result = 0;
      farthest_pptr = (std::max)(farthest_pptr, pptr());
      if (farthest_pptr && farthest_pptr > pbase()) {
        off_type delta = pptr() - farthest_pptr;
        int_type status = overflow();
        if (traits_type::eq_int_type(status, traits_type::eof())) result = -1;
        if (!py_seek.is_none()) py_seek(delta, 1);
      }
      else if (gptr() && gptr() < egptr()) {
        if (!py_seek.is_none()) py_seek(gptr() - egptr(), 1);
      }
      return result;
    }

    /// C.f. C++ standard section 27.5.2.4.2
    /** This implementation is optimised to look whether the position is within
        the buffers, so as to avoid calling Python seek or tell. It is
        important for many applications that the overhead of calling into Python
        is avoided as much as possible (e.g. parsers which may do a lot of
        backtracking)
    */
    virtual
    pos_type seekoff(off_type off, std::ios_base::seekdir way,
                     std::ios_base::openmode which=  std::ios_base::in
                                                   | std::ios_base::out)
    {
      /* In practice, "which" is either std::ios_base::in or out
         since we end up here because either seekp or seekg was called
         on the stream using this buffer. That simplifies the code
         in a few places.
      */
      int const failure = off_type(-1);

      if (py_seek.is_none()) {
        throw std::invalid_argument(
          "That Python file object has no 'seek' attribute");
      }

      // we need the read buffer to contain something!
      if (which == std::ios_base::in && !gptr()) {
        if (traits_type::eq_int_type(underflow(), traits_type::eof())) {
          return failure;
        }
      }

      // compute the whence parameter for Python seek
      int whence;
      switch (way) {
        case std::ios_base::beg:
          whence = 0;
          break;
        case std::ios_base::cur:
          whence = 1;
          break;
        case std::ios_base::end:
          whence = 2;
          break;
        default:
          return failure;
      }

      // Let's have a go
      off_type result;
      if (!seekoff_without_calling_python(off, way, which, result)) {
        // we need to call Python
        if (which == std::ios_base::out) overflow();
        if (way == std::ios_base::cur) {
          if      (which == std::ios_base::in)  off -= egptr() - gptr();
          else if (which == std::ios_base::out) off += pptr() - pbase();
        }
        py_seek(off, whence);
        result = off_type(py_tell().cast<off_type>());
        if (which == std::ios_base::in) underflow();
      }
      return result;
    }

    /// C.f. C++ standard section 27.5.2.4.2
    virtual
    pos_type seekpos(pos_type sp,
                     std::ios_base::openmode which=  std::ios_base::in
                                                   | std::ios_base::out)
    {
      return streambuf::seekoff(sp, std::ios_base::beg, which);
    }

  private:
    py::object py_read, py_write, py_seek, py_tell;

    std::size_t buffer_size;

    /* This is actually a Python bytes object and the actual read buffer is
       its internal data, i.e. an array of characters.
     */
    py::bytes read_buffer;

    /* A mere array of char's allocated on the heap at construction time and
       de-allocated only at destruction time.
    */
    char *write_buffer;

    off_type pos_of_read_buffer_end_in_py_file,
             pos_of_write_buffer_end_in_py_file;

    // the farthest place the buffer has been written into
    char *farthest_pptr;


    bool seekoff_without_calling_python(
      off_type off,
      std::ios_base::seekdir way,
      std::ios_base::openmode which,
      off_type & result)
    {
      // Buffer range and current position
      off_type buf_begin, buf_end, buf_cur, upper_bound;
      off_type pos_of_buffer_end_in_py_file;
      if (which == std::ios_base::in) {
        pos_of_buffer_end_in_py_file = pos_of_read_buffer_end_in_py_file;
        buf_begin = reinterpret_cast<std::streamsize>(eback());
        buf_cur = reinterpret_cast<std::streamsize>(gptr());
        buf_end = reinterpret_cast<std::streamsize>(egptr());
        upper_bound = buf_end;
      }
      else if (which == std::ios_base::out) {
        pos_of_buffer_end_in_py_file = pos_of_write_buffer_end_in_py_file;
        buf_begin = reinterpret_cast<std::streamsize>(pbase());
        buf_cur = reinterpret_cast<std::streamsize>(pptr());
        buf_end = reinterpret_cast<std::streamsize>(epptr());
        farthest_pptr = (std::max)(farthest_pptr, pptr());
        upper_bound = reinterpret_cast<std::streamsize>(farthest_pptr) + 1;
      }
      else {
           throw std::runtime_error(
             "Control flow passes through branch that should be unreachable.");
      }

      // Sought position in "buffer coordinate"
      off_type buf_sought;
      if (way == std::ios_base::cur) {
        buf_sought = buf_cur + off;
      }
      else if (way == std::ios_base::beg) {
        buf_sought = buf_end + (off - pos_of_buffer_end_in_py_file);
      }
      else if (way == std::ios_base::end) {
        return false;
      }
      else {
           throw std::runtime_error(
             "Control flow passes through branch that should be unreachable.");
      }

      // if the sought position is not in the buffer, give up
      if (buf_sought < buf_begin || buf_sought >= upper_bound) return false;

      // we are in wonderland
      if      (which == std::ios_base::in)  gbump(buf_sought - buf_cur);
      else if (which == std::ios_base::out) pbump(buf_sought - buf_cur);

      result = pos_of_buffer_end_in_py_file + (buf_sought - buf_end);
      return true;
    }

  public:

    class istream : public std::istream
    {
      public:
        istream(streambuf& buf) : std::istream(&buf)
        {
          exceptions(std::ios_base::badbit);
        }

        ~istream() { if (this->good()) this->sync(); }
    };

    class ostream : public std::ostream
    {
      public:
        ostream(streambuf& buf) : std::ostream(&buf)
        {
          exceptions(std::ios_base::badbit);
        }

        ~ostream() { if (this->good()) this->flush(); }
    };
};


struct streambuf_capsule
{
  streambuf python_streambuf;

  streambuf_capsule(
    py::object& python_file_obj,
    std::size_t buffer_size=0)
  :
    python_streambuf(python_file_obj, buffer_size)
  {}
};

struct ostream : private streambuf_capsule, streambuf::ostream
{
  ostream(
    py::object& python_file_obj,
    std::size_t buffer_size=0)
  :
    streambuf_capsule(python_file_obj, buffer_size),
    streambuf::ostream(python_streambuf)
  {}

  ~ostream()
  {
    if (this->good()){
      this->flush();
    }
  }
};

struct istream : private streambuf_capsule, streambuf::istream
{
  istream(
    py::object& python_file_obj,
    std::size_t buffer_size=0)
  :
    streambuf_capsule(python_file_obj, buffer_size),
    streambuf::istream(python_streambuf)
  {}

  ~istream()
  {
    if (this->good()) {
      this->sync();
    }
  }
};

};

namespace pybind11 { namespace detail {
    // Wrapping the stream objects with std::shared_ptr so that the streams can be used between API calls.
    // This enables one API call to open the stream and a follow-up call to read/write to it, with the stream
    // reference kept in C++.
    template <> struct type_caster<std::shared_ptr<pystream::istream>> {
    public:
        bool load(handle src, bool) {
            if (getattr(src, "read", py::none()).is_none()){
              return false;
            }

			obj = py::reinterpret_borrow<object>(src);
            value = std::shared_ptr<pystream::istream>(new pystream::istream(obj, 0));

            return true;
        }

    protected:
		py::object obj;
        std::shared_ptr<pystream::istream> value;

    public:
        static constexpr auto name = _("io.BytesIO");
        static handle cast(std::shared_ptr<pystream::istream> &src, return_value_policy policy, handle parent) {
            return none().release();
        }
        operator std::shared_ptr<pystream::istream>*() { return &value; }
        operator std::shared_ptr<pystream::istream>&() { return value; }
        template <typename _T> using cast_op_type = pybind11::detail::cast_op_type<_T>;
    };

    template <> struct type_caster<std::shared_ptr<pystream::ostream>> {
    public:
        bool load(handle src, bool) {
            if (getattr(src, "write", py::none()).is_none()){
              return false;
            }

			obj = py::reinterpret_borrow<object>(src);
            value = std::shared_ptr<pystream::ostream>(new pystream::ostream(obj, 0));

            return true;
        }

    protected:
		py::object obj;
        std::shared_ptr<pystream::ostream> value;

    public:
        static constexpr auto name = _("io.BytesIO");
        static handle cast(std::shared_ptr<pystream::ostream> &src, return_value_policy policy, handle parent) {
            return none().release();
        }
        operator std::shared_ptr<pystream::ostream>*() { return &value; }
        operator std::shared_ptr<pystream::ostream>&() { return value; }
        template <typename _T> using cast_op_type = pybind11::detail::cast_op_type<_T>;
    };
}} // namespace pybind11::detail

