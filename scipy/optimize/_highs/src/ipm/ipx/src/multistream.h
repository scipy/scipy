// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_MULTISTREAM_H_
#define IPX_MULTISTREAM_H_

#include <ostream>
#include <vector>

namespace ipx {

// Multistream allows redirecting output to multiple std::ostreams. It is used
// from Control (see control.h) for printing log messages simultaneously to
// standard output and a logfile.

class Multistream : public std::ostream {
public:
    Multistream() : std::ostream(nullptr) {
        std::ostream::rdbuf(&mbuffer_);
    }

    // adds a new stream to the object
    void add(std::ostream& os) {
        os.flush();
        mbuffer_.add(os.rdbuf());
    }

    // discards all streams
    void clear() {
        mbuffer_.clear();
    }

private:
    struct multibuffer : public std::streambuf {
        void add(std::streambuf* b) {
            buffers.push_back(b);
        }
        void clear() {
            buffers.clear();
        }
        int overflow(int c) override {
            for (std::streambuf* b : buffers)
                b->sputc(c);
            return c;
        }
    private:
        std::vector<std::streambuf*> buffers;
    };

    multibuffer mbuffer_;
};

}  // namespace ipx

#endif  // IPX_MULTISTREAM_H_
