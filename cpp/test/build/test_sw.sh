#!/bin/bash

cd allsky
cp ../*.nc .
./allsky rrtmgp-allsky.nc rrtmgp-data-sw-g224-2018-12-04.nc rrtmgp-cloud-optics-coeffs-sw.nc 1 1 >& tmpsw.txt || exit -1
UP=`grep flux_up tmpsw.txt | cut -d ":" -f 2`
DN=`grep flux_dn tmpsw.txt | grep -v flux_dn_dir | cut -d ":" -f 2`
DIR=`grep flux_dn_dir tmpsw.txt | cut -d ":" -f 2`

echo "flux_up     relative difference: `echo "($UP -19104.862129836212)/(19104.862129836212)*10^10" | bc -l` x 10^-10 "
echo "flux_dn     relative difference: `echo "($DN -38046.157649700355)/(38046.157649700355)*10^10" | bc -l` x 10^-10 "
echo "flux_dn_dir relative difference: `echo "($DIR-24998.593939345046)/(24998.593939345046)*10^10" | bc -l` x 10^-10 "

[[ "`echo "($UP -19104.862129836212)/(19104.862129836212) > 10^-10" | bc -l`" == "1" ]] && exit -1
[[ "`echo "($DN -38046.157649700355)/(38046.157649700355) > 10^-10" | bc -l`" == "1" ]] && exit -1
[[ "`echo "($DIR-24998.593939345046)/(24998.593939345046) > 10^-10" | bc -l`" == "1" ]] && exit -1
rm -f tmpsw.txt
 
exit 0

