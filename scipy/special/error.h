#pragma once


namespace scipy {
    namespace special {
#ifndef SP_SPECFUN_ERROR
        inline void set_error(const char *func_name, int code, const char *fmt, ...) {
            // nothing
        }
#else
        void set_error(const char *func_name, int code, const char *fmt, ...);
#endif
    }
}
