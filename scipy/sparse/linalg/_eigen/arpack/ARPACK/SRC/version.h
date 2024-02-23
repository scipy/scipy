/*

 In the current version, the parameter KAPPA in the Kahan's test
 for orthogonality is set to 0.717, the same as used by Gragg & Reichel.
 However computational experience indicates that this is a little too 
 strict and will frequently force reorthogonalization when it is not
 necessary to do so. 

 Also the "moving boundary" idea is not currently activated in the nonsymmetric
 code since it is not conclusive that it's the right thing to do all the time.  
 Requires further investigation.

 As of 02/01/93 Richard Lehoucq assumes software control of the codes from
 Phuong Vu. On 03/01/93 all the *.F files were migrated SCCS. The 1.1 version
 of codes are those received from Phuong Vu. The frozen version of 07/08/92
 is now considered version 1.1.

 Version 2.1 contains two new symmetric routines, sesrt and seupd. 
 Changes as well as bug fixes for version 1.1 codes that were only corrected 
 for programming bugs are version 1.2. These 1.2 versions will also be in version 2.1.
 Subroutine [d,s]saupd now requires slightly more workspace. See [d,s]saupd for the
 details. 

 \SCCS Information: @(#) 
  FILE: version.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2

 */

#define VERSION_NUMBER ' 2.1'
#define VERSION_DATE   ' 11/15/95'
