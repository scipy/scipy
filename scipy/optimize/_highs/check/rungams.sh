#!/bin/sh

# stop on error
set -e

if [ "$1" == "" ] ; then
  GAMSPATH=`dirname $(which gams)`
else
  GAMSPATH="$1"
fi

"$GAMSPATH/gamslib" trnsport

# try to print log if below there is an error, e.g., gams fails or grep does not find
trap "cat trnsport.log" ERR

"$GAMSPATH/gams" trnsport.gms LP=HIGHS LO=2

grep "**** SOLVER STATUS     1 Normal Completion" trnsport.lst
grep "**** MODEL STATUS      1 Optimal" trnsport.lst
grep "**** OBJECTIVE VALUE              153.6750" trnsport.lst

rm -f trnsport.{gms,log,lst}
