#! /bin/csh

set ofile = install.out			# output file

echo '---- SINGLE PRECISION' >! $ofile
./testslamch >> $ofile
echo '' >> $ofile
echo ---- DOUBLE PRECISION >> $ofile
./testdlamch >> $ofile
echo '' >> $ofile
echo ---- TIMER >> $ofile
./testtimer >> $ofile


