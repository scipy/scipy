//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Matei David (matei@cs.toronto.edu)
//---------------------------------------------------------

// Reference:
// http://stackoverflow.com/questions/14086417/how-to-write-custom-input-stream-in-c

#pragma once

#include <zlib.h>

#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#include "../extern/zstr/strict_fstream.hpp"

namespace zstr {

static const std::size_t default_buff_size = static_cast<std::size_t>(1 << 20);

/// Exception class thrown by failed zlib operations.
class Exception : public std::ios_base::failure {
 public:
  static std::string error_to_message(z_stream* zstrm_p, int ret) {
    std::string msg = "zlib: ";
    switch (ret) {
      case Z_STREAM_ERROR:
        msg += "Z_STREAM_ERROR: ";
        break;
      case Z_DATA_ERROR:
        msg += "Z_DATA_ERROR: ";
        break;
      case Z_MEM_ERROR:
        msg += "Z_MEM_ERROR: ";
        break;
      case Z_VERSION_ERROR:
        msg += "Z_VERSION_ERROR: ";
        break;
      case Z_BUF_ERROR:
        msg += "Z_BUF_ERROR: ";
        break;
      default:
        std::ostringstream oss;
        oss << ret;
        msg += "[" + oss.str() + "]: ";
        break;
    }
    if (zstrm_p->msg) {
      msg += zstrm_p->msg;
    }
    msg +=
        " ("
        "next_in: " +
        std::to_string(uintptr_t(zstrm_p->next_in)) +
        ", avail_in: " + std::to_string(uintptr_t(zstrm_p->avail_in)) +
        ", next_out: " + std::to_string(uintptr_t(zstrm_p->next_out)) +
        ", avail_out: " + std::to_string(uintptr_t(zstrm_p->avail_out)) + ")";
    return msg;
  }

  Exception(z_stream* zstrm_p, int ret)
      : std::ios_base::failure(error_to_message(zstrm_p, ret)) {}
};  // class Exception

namespace detail {

class z_stream_wrapper : public z_stream {
 public:
  z_stream_wrapper(bool _is_input, int _level, int _window_bits)
      : is_input(_is_input) {
    this->zalloc = nullptr;  // Z_NULL
    this->zfree = nullptr;   // Z_NULL
    this->opaque = nullptr;  // Z_NULL
    int ret;
    if (is_input) {
      this->avail_in = 0;
      this->next_in = nullptr;  // Z_NULL
      ret = inflateInit2(this, _window_bits ? _window_bits : 15 + 32);
    } else {
      ret = deflateInit2(this, _level, Z_DEFLATED,
                         _window_bits ? _window_bits : 15 + 16, 8,
                         Z_DEFAULT_STRATEGY);
    }
    if (ret != Z_OK) throw Exception(this, ret);
  }
  ~z_stream_wrapper() {
    if (is_input) {
      inflateEnd(this);
    } else {
      deflateEnd(this);
    }
  }

 private:
  bool is_input;
};  // class z_stream_wrapper

}  // namespace detail

class istreambuf : public std::streambuf {
 public:
  istreambuf(std::streambuf* _sbuf_p,
             std::size_t _buff_size = default_buff_size,
             bool _auto_detect = true, int _window_bits = 0)
      : sbuf_p(_sbuf_p),
        in_buff(),
        in_buff_start(nullptr),
        in_buff_end(nullptr),
        out_buff(),
        zstrm_p(nullptr),
        buff_size(_buff_size),
        auto_detect(_auto_detect),
        auto_detect_run(false),
        is_text(false),
        window_bits(_window_bits) {
    assert(sbuf_p);
    in_buff = std::unique_ptr<char[]>(new char[buff_size]);
    in_buff_start = in_buff.get();
    in_buff_end = in_buff.get();
    out_buff = std::unique_ptr<char[]>(new char[buff_size]);
    setg(out_buff.get(), out_buff.get(), out_buff.get());
  }

  istreambuf(const istreambuf&) = delete;
  istreambuf& operator=(const istreambuf&) = delete;

  pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                   std::ios_base::openmode which) override {
    if (off != 0 || dir != std::ios_base::cur) {
      return std::streambuf::seekoff(off, dir, which);
    }

    if (!zstrm_p) {
      return 0;
    }

    return static_cast<long int>(zstrm_p->total_out -
                                 static_cast<uLong>(in_avail()));
  }

  std::streambuf::int_type underflow() override {
    if (this->gptr() == this->egptr()) {
      // pointers for free region in output buffer
      char* out_buff_free_start = out_buff.get();
      int tries = 0;
      do {
        if (++tries > 1000) {
          throw std::ios_base::failure(
              "Failed to fill buffer after 1000 tries");
        }

        // read more input if none available
        if (in_buff_start == in_buff_end) {
          // empty input buffer: refill from the start
          in_buff_start = in_buff.get();
          std::streamsize sz = sbuf_p->sgetn(
              in_buff.get(), static_cast<std::streamsize>(buff_size));
          in_buff_end = in_buff_start + sz;
          if (in_buff_end == in_buff_start) break;  // end of input
        }
        // auto detect if the stream contains text or deflate data
        if (auto_detect && !auto_detect_run) {
          auto_detect_run = true;
          unsigned char b0 = *reinterpret_cast<unsigned char*>(in_buff_start);
          unsigned char b1 =
              *reinterpret_cast<unsigned char*>(in_buff_start + 1);
          // Ref:
          // http://en.wikipedia.org/wiki/Gzip
          // http://stackoverflow.com/questions/9050260/what-does-a-zlib-header-look-like
          is_text = !(in_buff_start + 2 <= in_buff_end &&
                      ((b0 == 0x1F && b1 == 0x8B)     // gzip header
                       || (b0 == 0x78 && (b1 == 0x01  // zlib header
                                          || b1 == 0x9C || b1 == 0xDA))));
        }
        if (is_text) {
          // simply swap in_buff and out_buff, and adjust pointers
          assert(in_buff_start == in_buff.get());
          std::swap(in_buff, out_buff);
          out_buff_free_start = in_buff_end;
          in_buff_start = in_buff.get();
          in_buff_end = in_buff.get();
        } else {
          // run inflate() on input
          if (!zstrm_p)
            zstrm_p = std::unique_ptr<detail::z_stream_wrapper>(
                new detail::z_stream_wrapper(true, Z_DEFAULT_COMPRESSION,
                                             window_bits));
          zstrm_p->next_in =
              reinterpret_cast<decltype(zstrm_p->next_in)>(in_buff_start);
          zstrm_p->avail_in = uint32_t(in_buff_end - in_buff_start);
          zstrm_p->next_out = reinterpret_cast<decltype(zstrm_p->next_out)>(
              out_buff_free_start);
          zstrm_p->avail_out =
              uint32_t((out_buff.get() + buff_size) - out_buff_free_start);
          int ret = inflate(zstrm_p.get(), Z_NO_FLUSH);
          // process return code
          if (ret != Z_OK && ret != Z_STREAM_END)
            throw Exception(zstrm_p.get(), ret);
          // update in&out pointers following inflate()
          in_buff_start =
              reinterpret_cast<decltype(in_buff_start)>(zstrm_p->next_in);
          in_buff_end = in_buff_start + zstrm_p->avail_in;
          out_buff_free_start = reinterpret_cast<decltype(out_buff_free_start)>(
              zstrm_p->next_out);
          assert(out_buff_free_start + zstrm_p->avail_out ==
                 out_buff.get() + buff_size);

          if (ret == Z_STREAM_END) {
            // if stream ended, deallocate inflator
            zstrm_p.reset();
          }
        }
      } while (out_buff_free_start == out_buff.get());
      // 2 exit conditions:
      // - end of input: there might or might not be output available
      // - out_buff_free_start != out_buff: output available
      this->setg(out_buff.get(), out_buff.get(), out_buff_free_start);
    }
    return this->gptr() == this->egptr()
               ? traits_type::eof()
               : traits_type::to_int_type(*this->gptr());
  }

 private:
  std::streambuf* sbuf_p;
  std::unique_ptr<char[]> in_buff;
  char* in_buff_start;
  char* in_buff_end;
  std::unique_ptr<char[]> out_buff;
  std::unique_ptr<detail::z_stream_wrapper> zstrm_p;
  std::size_t buff_size;
  bool auto_detect;
  bool auto_detect_run;
  bool is_text;
  int window_bits;

};  // class istreambuf

class ostreambuf : public std::streambuf {
 public:
  ostreambuf(std::streambuf* _sbuf_p,
             std::size_t _buff_size = default_buff_size,
             int _level = Z_DEFAULT_COMPRESSION, int _window_bits = 0)
      : sbuf_p(_sbuf_p),
        in_buff(),
        out_buff(),
        zstrm_p(new detail::z_stream_wrapper(false, _level, _window_bits)),
        buff_size(_buff_size) {
    assert(sbuf_p);
    in_buff = std::unique_ptr<char[]>(new char[buff_size]);
    out_buff = std::unique_ptr<char[]>(new char[buff_size]);
    setp(in_buff.get(), in_buff.get() + buff_size);
  }

  ostreambuf(const ostreambuf&) = delete;
  ostreambuf& operator=(const ostreambuf&) = delete;

  int deflate_loop(int flush) {
    while (true) {
      zstrm_p->next_out =
          reinterpret_cast<decltype(zstrm_p->next_out)>(out_buff.get());
      zstrm_p->avail_out = uint32_t(buff_size);
      int ret = deflate(zstrm_p.get(), flush);
      if (ret != Z_OK && ret != Z_STREAM_END && ret != Z_BUF_ERROR) {
        failed = true;
        throw Exception(zstrm_p.get(), ret);
      }
      std::streamsize sz = sbuf_p->sputn(
          out_buff.get(),
          reinterpret_cast<decltype(out_buff.get())>(zstrm_p->next_out) -
              out_buff.get());
      if (sz != reinterpret_cast<decltype(out_buff.get())>(zstrm_p->next_out) -
                    out_buff.get()) {
        // there was an error in the sink stream
        return -1;
      }
      if (ret == Z_STREAM_END || ret == Z_BUF_ERROR || sz == 0) {
        break;
      }
    }
    return 0;
  }

  virtual ~ostreambuf() {
    // flush the zlib stream
    //
    // NOTE: Errors here (sync() return value not 0) are ignored, because we
    // cannot throw in a destructor. This mirrors the behaviour of
    // std::basic_filebuf::~basic_filebuf(). To see an exception on error,
    // close the ofstream with an explicit call to close(), and do not rely
    // on the implicit call in the destructor.
    //
    if (!failed) try {
        sync();
      } catch (...) {
      }
  }
  std::streambuf::int_type overflow(
      std::streambuf::int_type c = traits_type::eof()) override {
    zstrm_p->next_in = reinterpret_cast<decltype(zstrm_p->next_in)>(pbase());
    zstrm_p->avail_in = uint32_t(pptr() - pbase());
    while (zstrm_p->avail_in > 0) {
      int r = deflate_loop(Z_NO_FLUSH);
      if (r != 0) {
        setp(nullptr, nullptr);
        return traits_type::eof();
      }
    }
    setp(in_buff.get(), in_buff.get() + buff_size);
    return traits_type::eq_int_type(c, traits_type::eof())
               ? traits_type::eof()
               : sputc(char_type(c));
  }
  int sync() override {
    // first, call overflow to clear in_buff
    overflow();
    if (!pptr()) return -1;
    // then, call deflate asking to finish the zlib stream
    zstrm_p->next_in = nullptr;
    zstrm_p->avail_in = 0;
    if (deflate_loop(Z_FINISH) != 0) return -1;
    deflateReset(zstrm_p.get());
    return 0;
  }

 private:
  std::streambuf* sbuf_p = nullptr;
  std::unique_ptr<char[]> in_buff;
  std::unique_ptr<char[]> out_buff;
  std::unique_ptr<detail::z_stream_wrapper> zstrm_p;
  std::size_t buff_size;
  bool failed = false;

};  // class ostreambuf

class istream : public std::istream {
 public:
  istream(std::istream& is, std::size_t _buff_size = default_buff_size,
          bool _auto_detect = true, int _window_bits = 0)
      : std::istream(new istreambuf(is.rdbuf(), _buff_size, _auto_detect,
                                    _window_bits)) {
    exceptions(std::ios_base::badbit);
  }
  explicit istream(std::streambuf* sbuf_p)
      : std::istream(new istreambuf(sbuf_p)) {
    exceptions(std::ios_base::badbit);
  }
  virtual ~istream() { delete rdbuf(); }
};  // class istream

class ostream : public std::ostream {
 public:
  ostream(std::ostream& os, std::size_t _buff_size = default_buff_size,
          int _level = Z_DEFAULT_COMPRESSION, int _window_bits = 0)
      : std::ostream(
            new ostreambuf(os.rdbuf(), _buff_size, _level, _window_bits)) {
    exceptions(std::ios_base::badbit);
  }
  explicit ostream(std::streambuf* sbuf_p)
      : std::ostream(new ostreambuf(sbuf_p)) {
    exceptions(std::ios_base::badbit);
  }
  virtual ~ostream() { delete rdbuf(); }
};  // class ostream

namespace detail {

template <typename FStream_Type>
struct strict_fstream_holder {
  strict_fstream_holder(const std::string& filename,
                        std::ios_base::openmode mode = std::ios_base::in)
      : _fs(filename, mode) {}
  strict_fstream_holder() = default;
  FStream_Type _fs{};
};  // class strict_fstream_holder

}  // namespace detail

class ifstream
    : private detail::strict_fstream_holder<strict_fstream::ifstream>,
      public std::istream {
 public:
  explicit ifstream(const std::string filename,
                    std::ios_base::openmode mode = std::ios_base::in,
                    size_t buff_size = default_buff_size)
      : detail::strict_fstream_holder<strict_fstream::ifstream>(
            filename, mode
#ifdef _WIN32  // to avoid problems with conversion of \r\n, only windows as
               // otherwise there are problems on mac
                          | std::ios_base::binary
#endif
            ),
        std::istream(new istreambuf(_fs.rdbuf(), buff_size)) {
    exceptions(std::ios_base::badbit);
  }
  explicit ifstream()
      : detail::strict_fstream_holder<strict_fstream::ifstream>(),
        std::istream(new istreambuf(_fs.rdbuf())) {}
  void close() { _fs.close(); }
  void open(const std::string filename,
            std::ios_base::openmode mode = std::ios_base::in) {
    _fs.open(filename, mode
#ifdef _WIN32  // to avoid problems with conversion of \r\n, only windows as
               // otherwise there are problems on mac
                           | std::ios_base::binary
#endif
    );
    // make sure the previous buffer is deleted by putting it into a unique_ptr
    // and set a new one after opening file
    std::unique_ptr<std::streambuf> oldbuf(rdbuf(new istreambuf(_fs.rdbuf())));
    // call move assignment operator on istream which does not alter the stream
    // buffer
    std::istream::operator=(std::istream(rdbuf()));
  }
  bool is_open() const { return _fs.is_open(); }
  virtual ~ifstream() {
    if (_fs.is_open()) close();
    if (rdbuf()) delete rdbuf();
  }

  /// Return the position within the compressed file (wrapped filestream)
  std::streampos compressed_tellg() { return _fs.tellg(); }
};  // class ifstream

class ofstream
    : private detail::strict_fstream_holder<strict_fstream::ofstream>,
      public std::ostream {
 public:
  explicit ofstream(const std::string filename,
                    std::ios_base::openmode mode = std::ios_base::out,
                    int level = Z_DEFAULT_COMPRESSION,
                    size_t buff_size = default_buff_size)
      : detail::strict_fstream_holder<strict_fstream::ofstream>(
            filename, mode | std::ios_base::binary),
        std::ostream(new ostreambuf(_fs.rdbuf(), buff_size, level)) {
    exceptions(std::ios_base::badbit);
  }
  explicit ofstream()
      : detail::strict_fstream_holder<strict_fstream::ofstream>(),
        std::ostream(new ostreambuf(_fs.rdbuf())) {}
  void close() {
    std::ostream::flush();
    _fs.close();
  }
  void open(const std::string filename,
            std::ios_base::openmode mode = std::ios_base::out,
            int level = Z_DEFAULT_COMPRESSION) {
    flush();
    _fs.open(filename, mode | std::ios_base::binary);
    std::ostream::operator=(
        std::ostream(new ostreambuf(_fs.rdbuf(), default_buff_size, level)));
  }
  bool is_open() const { return _fs.is_open(); }
  ofstream& flush() {
    std::ostream::flush();
    _fs.flush();
    return *this;
  }
  virtual ~ofstream() {
    if (_fs.is_open()) close();
    if (rdbuf()) delete rdbuf();
  }

  // Return the position within the compressed file (wrapped filestream)
  std::streampos compressed_tellp() { return _fs.tellp(); }
};  // class ofstream

}  // namespace zstr
