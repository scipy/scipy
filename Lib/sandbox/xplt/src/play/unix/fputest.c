/*
 * fputest.c -- $Id$
 * test for SIGFPE delivery
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#ifndef NO_FPUSET
#define FROM_FPUTEST
#include "fpuset.c"
#else
extern void u_fpu_setup(int when);
void u_fpu_setup(int when) {}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include <setjmp.h>
static jmp_buf u_jmp_target;

extern void u_sigfpe(int sig);
extern double reciprocal(double x);
extern double quotient(double x, double y);
extern double powpow(double x, int n);
extern int decrement(int *ip);
static int always_true = 0;

int
main(int argc, char *argv[])
{
  int i = 5;
  double zero = (double)(argc>1000);
  double huge = (argc>1000)? 2.0 : 1.e20;
  always_true = !(argc>1000);

  u_fpu_setup(-1);
  if (setjmp(u_jmp_target)) u_fpu_setup(0);

  /* signal *ought* to be enough to get SIGFPE delivered
   * -- but it never is -- see README.fpu */
  signal(SIGFPE, &u_sigfpe);

  /* need to make sure that loop index i actually decrements
   * despite interrupt */
  while (decrement(&i)) {
    printf("SIGFPE handling failed on pass %d, 1./0. = %g\n", 5-i,
           reciprocal(zero));
    if (i) break;
  }
  if (!i) {
    i = 5;
    if (setjmp(u_jmp_target)) u_fpu_setup(0);
    while (decrement(&i)) {
      printf("SIGFPE handling failed on pass %d, 0./0. = %g\n", 5-i,
             quotient(zero, zero));
      if (i) break;
    }
    if (!i) {
      i = 5;
      if (setjmp(u_jmp_target)) u_fpu_setup(0);
      while (decrement(&i)) {
        printf("SIGFPE handling failed on pass %d, 10.^20480 = %g\n", 5-i,
               powpow(huge, 10));
        if (i) break;
      }
      if (!i) {
        if (setjmp(u_jmp_target)) {
          puts("SIGFPE improperly generated on underflow");
          i = 11;
        } else {
          double x = powpow(1./huge, 10);
          if (x != 0.0)
            printf("SIGFPE handling works, but 10.^-20480 = %g\n", x);
          else
            puts("SIGFPE handling works properly");
        }
      }
    }
  }
  exit(i? 2 : 0);
}

void
u_sigfpe(int sig)
{
  if (sig==SIGFPE) {
    u_fpu_setup(1);
    signal(SIGFPE, &u_sigfpe);
    longjmp(u_jmp_target, 1);
  } else {
    puts("u_sigfpe called, but with bad parameter");
  }
  exit(1);
}

int
decrement(int *ip)
{
  int i = *ip;
  if (always_true) i--;
  else i -= 2;
  *ip = i;
  return i>0;
}

double
reciprocal(double x)
{
  return 1./x;
}

double
quotient(double x, double y)
{
  return x/y;
}

double
powpow(double x, int n)
{
  double y = x;
  while (n--) y *= y;
  return y;
}
