/*
 * mfcmain.cpp -- $Id$
 * MFC main program stub
 *
 * Copyright (c) 2000.  See accompanying LEGAL file for details.
 */

extern "C" {
  extern int on_launch(int argc, char *argv[]);
}
#include "mfcapp.h"

mfc_boss the_boss(on_launch);
