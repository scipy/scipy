/*
 * fpuset.c -- $Id$
 * set up FPU to trap floating point exceptions
 * - this is very non-portable, not covered by ANSI C, POSIX, or even C9X
 * - if you port to a new platform (eg- Ultrix) please contact the author
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#ifndef FROM_FPUTEST
# include "config.h"
# include "playu.h"
#else
extern void u_fpu_setup(int when);
#endif

/* when = -1 for initial call before setjmp
 *         0 after each longjmp out of interrupt handler (after setjmp)
 *         1 inside interrupt handler before signal() re-enables SIGFPE
 */

#if defined(FPU_DIGITAL) || defined(FPU_ALPHA_LINUX)

/* FPU_ALPHA_LINUX: see README.fpu */
/* man pages: exception_intro, ieee */
# ifdef FPU_DIGITAL
#  include <machine/fpu.h>
# else
   extern void ieee_set_fp_control(long);
#  define IEEE_TRAP_ENABLE_INV 0x000002
#  define IEEE_TRAP_ENABLE_DZE 0x000004
#  define IEEE_TRAP_ENABLE_OVF 0x000008
#  define IEEE_MAP_DMZ         (1UL<<12)
#  define IEEE_MAP_UMZ         (1UL<<13)
# endif
void
u_fpu_setup(int when)
{
  /* possibly should include IEEE_MAP_DMZ and IEEE_MAP_UMZ
   * to map denorm inputs and underflowed outputs to zero
   * --however, these apparently only have an effect for software
   * completed operations (the hardware always maps underflows to zero)
   */
  if (when < 0) {
    ieee_set_fp_control(IEEE_TRAP_ENABLE_INV | IEEE_TRAP_ENABLE_DZE |
                        IEEE_TRAP_ENABLE_OVF);
  }
}

#elif defined(FPU_AIX)

/* man pages: fp_trap, fp_enable */
#include <fptrap.h>
void
u_fpu_setup(int when)
{
  if (when) {
    fp_trap(FP_TRAP_FASTMODE);
    fp_enable(TRP_INVALID | TRP_DIV_BY_ZERO | TRP_OVERFLOW);
  }
}

#elif defined(FPU_HPUX)

/* man pages: fpsetmask
 * library: -lm */
/* HPUX turns off FP_X_* without this (_INCLUDE_HPUX_SOURCE) */
#ifndef _HPUX_SOURCE
#define _HPUX_SOURCE 1
#endif
#include <math.h>
void
u_fpu_setup(int when)
{
  if (when <= 0) {
    fpsetmask(FP_X_INV | FP_X_DZ | FP_X_OFL);  /* 0x1c */
    fpsetfastmode(1);    /* fast underflows */
  }
}

#elif defined(FPU_IRIX)

/* man pages: handle_sigfpes, note lethal TRAP_FPE environment variable
 * library: -lfpe
 * note: earlier versions used get_fpc_csr/set_fpc_csr?, sys/fpu.h */
#include <sigfpe.h>
void
u_fpu_setup(int when)
{
  if (when < 0) {
    extern void u_sigfpe(int sig);  /* from handler.c (or fputest.c) */
    handle_sigfpes(_ON, _EN_OVERFL|_EN_DIVZERO|_EN_INVALID,
                   (void (*)())0, _USER_HANDLER, (void (*)())&u_sigfpe);
  }
}

#elif defined(FPU_SOLARIS)

/* man pages: fpsetmask
 *    Sun's -fnonstd compiler switch switches between __fnonstd.o
 *      and __fstd.o under Solaris, as far as I can tell.  Use FPU_IGNORE
 *        if you do this.  */
#include <ieeefp.h>
void
u_fpu_setup(int when)
{
  if (when < 0) {
    fpsetmask(FP_X_INV | FP_X_DZ | FP_X_OFL);
    /* this doesn't set the "nonstandard arithmetic" bit, which prevents
     * software emulation of IEEE gradual underflow
     * -- apparently no way to do this in libc (see FPU_GCC_SPARC) */
  }
}

#elif defined(FPU_SUN4)

/* man pages: ieee_handler
 *               nonstandard_arithmetic is undocumented, but rumored
 *               to be important to get rapid underflows
 * library: -lsunmath (under /usr/lang hierarchy)
 *          may also be in -lm (standard libm)?
 *   note: libsunmath.a is provided by Sun only if you purchase their
 *         compilers; if you are trying to compile with gcc on a SPARC
 *         architecture, try FPU_GCC_SPARC
 *    Sun's -fnonstd compiler switch buggers crt1.o under SunOS 4,
 *      as far as I can tell.  Use FPU_IGNORE if you do this
 *      (not possible with gcc?).  */
void
u_fpu_setup(int when)
{
  if (when < 0) {
    extern void u_sigfpe(int sig);  /* from handler.c (or fputest.c) */
    nonstandard_arithmetic();
    ieee_handler("set","common", &u_sigfpe);
  }
}

#elif defined(FPU_UNICOS)

/* delivers SIGFPE by default, this just arranges to trap on
 * libm errors as well */
void
u_fpu_setup(int when)
{
  if (when < 0) {
    int flag = -1;
    libmset(&flag);
  }
}

#elif defined(FPU_GNU_FENV)

/* GCC enhanced C9X fenv.h interface by adding feenableexcept */
#include <fenv.h>
void
u_fpu_setup(int when)
{
  if (when <= 0)
    feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
}

#elif defined(FPU_GCC_I86)

/* see also: fpu_control.h or i386/fpu_control.h, __setfpucw function */
void
u_fpu_setup(int when)
{
  if (when) {
    unsigned int fpucw = 0x1372;
    __asm__ ("fldcw %0" : : "m" (fpucw));
  }
}

#elif defined(FPU_GCC_POWERPC)

void
u_fpu_setup(int when)
{
  if (when) {
    unsigned int tmp[2] __attribute__ ((__aligned__(8)));
    tmp[0] = 0xFFF80000; /* More-or-less arbitrary; this is a QNaN. */
    tmp[1] = 0xd0;
    __asm__ ("lfd 0,%0; mtfsf 255,0" : : "m" (*tmp) : "fr0");
  }
}

#elif defined(FPU_GCC_SPARC)

void
u_fpu_setup(int when)
{
  if (when < 0) {
    unsigned int fpucw = 0xd400000;  /* the 4 is nonstandard arithmetic bit */
    __asm__ ("ld %0,%%fsr" : : "m" (fpucw));
  }
}

#elif defined(FPU_GCC_M68K)

/* works on NeXT as well as m68k Linux */
void
u_fpu_setup(int when)
{
  if (when <= 0) {
    asm("fmovel     #0x7400,fpcr");   /* set OVFL and ZD bits */
    /* unsigned int fpucw = 0x7400;
     * __asm__ volatile ("fmove%.l %0, %!" : : "dm" (fpucw)); */
    /* includes bit to trap on signalling NaN (may affect libm behavior) */
  }
}

#elif defined(FPU_GCC_ARM)

void
u_fpu_setup(int when)
{
  if (when <= 0) {
    unsigned int fpucw = 0x70200;
    __asm__ ("wfs %0" : : "r" (fpucw));
    /* includes bit to trap on signalling NaN (may affect libm behavior) */
  }
}

#elif defined(FPU_IGNORE)

void
u_fpu_setup(int when)
{
}

#elif defined(FPU_MACOSX)

#include <architecture/ppc/fp_regs.h>
#include <mach/mach.h>
#include <pthread.h>

static void *fpu_fpe_enable(void *arg);
#define FE0_MASK (1<<11)
#define FE1_MASK (1<<8)
/* FE0  FE1   exceptions enabled if either FE0 or FE1 set
 *  0    0    -- floating-point exceptions disabled
 *  0    1    -- floating-point imprecise nonrecoverable
 *  1    0    -- floating-point imprecise recoverable
 *  1    1    -- floating-point precise mode
 */
/* for Darwin version 6.0 (MacOS X 10.2) FE0=FE1=1 initially
 * for Darwin version 5.5 (MacOS X <=10.1) FE0=FE1=0 initially
 * Darwin 5.5 resets MSR to FE0=FE1=0 after each SIGFPE
 * Darwin 6.0 does not reset MSR?  leave MSR reset code in case?
 */

/* a thread cannot get or set its own MSR bits */
static void *
fpu_fpe_enable(void *arg)
{
  thread_t t = *(thread_t *)arg;
  struct ppc_thread_state state;
  unsigned int state_size = PPC_THREAD_STATE_COUNT;
  if (thread_get_state(t, PPC_THREAD_STATE,
                       (natural_t *)&state, &state_size) == KERN_SUCCESS) {
    state.srr1 |= FE1_MASK;
    state.srr1 &= ~FE0_MASK;
    thread_set_state(t, PPC_THREAD_STATE, (natural_t *)&state, state_size);
  }
  return 0;
}

void
u_fpu_setup(int when)
{
  static volatile int looping = 0;
  if (when) {
    ppc_fp_scr_t r = get_fp_scr();
    /* turn off exception bits to prevent immediate re-fault */
    r.fx = r.fex = r.vx = r.ox = r.ux = r.zx = r.xx = r.vx_snan = r.vx_isi =
      r.vx_idi = r.vx_zdz = r.vx_imz = r.vx_xvc = r.vx_cvi = r.vx_soft = 0;
    /* these only have to be set once, but may as well set anyway */
    r.ve = 1;  /* invalid */
    r.oe = 1;  /* overflow */
    r.ue = 0;  /* underflow */
    r.ze = 1;  /* zero divide */
    r.xe = 0;  /* inexact */
    if (!looping) {
      looping |= 1;
      set_fp_scr(r);
      looping &= ~1;
    }
  }
  if (when <= 0) {
    thread_t self = mach_thread_self();
    pthread_t enabler;
    if (!looping) {
      looping |= 2;
      if (!pthread_create(&enabler, 0, fpu_fpe_enable, &self))
        pthread_join(enabler, 0);
      looping &= ~2;
    }
  }
  looping = 0;
}

#else

#error <read play/unix/README.fpu for help>

#endif
