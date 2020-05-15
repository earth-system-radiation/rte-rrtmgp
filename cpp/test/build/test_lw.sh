#!/bin/bash

cd allsky
cp ../*.nc .
./allsky rrtmgp-allsky.nc rrtmgp-data-lw-g256-2018-12-04.nc rrtmgp-cloud-optics-coeffs-lw.nc 1 1 >& tmplw.txt || exit -1
UP=`grep flux_up tmplw.txt | cut -d ":" -f 2`
DN=`grep flux_dn tmplw.txt | cut -d ":" -f 2`

echo "flux_up relative difference: `echo "($UP-10264.518998579415)/(10264.518998579415)*10^10" | bc -l` x 10^-10 "
echo "flux_dn relative difference: `echo "($DN-6853.2350138542843)/(6853.2350138542843)*10^10" | bc -l` x 10^-10 "

[[ "`echo "($UP-10264.518998579415)/(10264.518998579415) > 10^-10" | bc -l`" == "1" ]] && exit -1
[[ "`echo "($DN-6853.2350138542843)/(6853.2350138542843) > 10^-10" | bc -l`" == "1" ]] && exit -1
rm -f tmplw.txt

exit 0

