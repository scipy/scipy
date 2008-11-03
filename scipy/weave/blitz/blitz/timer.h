/***************************************************************************
 * blitz/Timer.h        Timer class, for use in benchmarking
 *
 * $Id: timer.h 1414 2005-11-01 22:04:59Z cookedm $
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This code was relicensed under the modified BSD license for use in SciPy
 * by Todd Veldhuizen (see LICENSE.txt in the weave directory).
 *
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ***************************************************************************/

// This class is not portable to non System V platforms.
// It will need to be rewritten for Windows, NT, Mac.
// NEEDS_WORK

#ifndef BZ_TIMER_H
#define BZ_TIMER_H

#ifndef BZ_BLITZ_H
 #include <blitz/blitz.h>
#endif

#ifdef BZ_HAVE_RUSAGE
 #include <sys/resource.h>
#else
 #include <time.h>
#endif

BZ_NAMESPACE(blitz)

class Timer {

public:
    Timer() 
    { 
        state_ = uninitialized;
    }

    void start()
    { 
        state_ = running;
        t1_ = systemTime();
    }

    void stop()
    {
        t2_ = systemTime();
        BZPRECONDITION(state_ == running);
        state_ = stopped;
    }
    
/* Compaq cxx compiler in ansi mode cannot print out long double type! */
#if defined(__DECCXX)
    double elapsedSeconds()
#else
    long double elapsedSeconds()
#endif
    {
        BZPRECONDITION(state_ == stopped);
        return t2_ - t1_;
    }

private:
    Timer(Timer&) { }
    void operator=(Timer&) { }

    long double systemTime()
    {
#ifdef BZ_HAVE_RUSAGE
        getrusage(RUSAGE_SELF, &resourceUsage_);
        double seconds = resourceUsage_.ru_utime.tv_sec 
            + resourceUsage_.ru_stime.tv_sec;
        double micros  = resourceUsage_.ru_utime.tv_usec 
            + resourceUsage_.ru_stime.tv_usec;
        return seconds + micros/1.0e6;
#else
        return clock() / (long double) CLOCKS_PER_SEC;
#endif
    }

    enum { uninitialized, running, stopped } state_;

#ifdef BZ_HAVE_RUSAGE
    struct rusage resourceUsage_;
#endif

    long double t1_, t2_;
};

BZ_NAMESPACE_END

#endif // BZ_TIMER_H

