
#ifndef LIBHIGHS_EXPORT_H
#define LIBHIGHS_EXPORT_H

#ifdef LIBHIGHS_STATIC_DEFINE
#  define LIBHIGHS_EXPORT
#  define LIBHIGHS_NO_EXPORT
#else
#  ifndef LIBHIGHS_EXPORT
#    ifdef libhighs_EXPORTS
        /* We are building this library */
#        if defined(_MSC_VER)
#          define LIBHIGHS_EXPORT __declspec(dllexport)
#        else
#          define LIBHIGHS_EXPORT __attribute__((visibility("default")))
#        endif
#    else
        /* We are using this library */
#        if defined(_MSC_VER)
#          define LIBHIGHS_EXPORT __declspec(dllexport)
#        else
#          define LIBHIGHS_EXPORT __attribute__((visibility("default")))
#        endif
#    endif
#  endif

#  ifndef LIBHIGHS_NO_EXPORT
#    define LIBHIGHS_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef LIBHIGHS_DEPRECATED
#  define LIBHIGHS_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef LIBHIGHS_DEPRECATED_EXPORT
#  define LIBHIGHS_DEPRECATED_EXPORT LIBHIGHS_EXPORT LIBHIGHS_DEPRECATED
#endif

#ifndef LIBHIGHS_DEPRECATED_NO_EXPORT
#  define LIBHIGHS_DEPRECATED_NO_EXPORT LIBHIGHS_NO_EXPORT LIBHIGHS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef LIBHIGHS_NO_DEPRECATED
#    define LIBHIGHS_NO_DEPRECATED
#  endif
#endif

#endif /* LIBHIGHS_EXPORT_H */
