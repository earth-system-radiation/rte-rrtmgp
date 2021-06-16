. ~/python-venv/elastic/bin/activate

module load intel/17.0.6 
export FC="ifort"

export FCFLAGS="-m64 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132 -check bounds,uninit,pointers,stack -stand f08"
(NO OPTIMIZATION FLAGS)

export FCINCLUDE="-I/work/k20200/k202144/ICON/Reordering/Develop_brange/rte-rrtmgp"

export RRTMGP_ROOT="/work/k20200/k202144/ICON/Reordering/Develop_brange/rte-rrtmgp"

export NFHOME="/sw/spack-rhel6/netcdf-fortran-4.5.3-nqwtm5"
export NCHOME="/sw/spack-rhel6/netcdf-c-4.7.4-hlp365"

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NFHOME/lib:$NCHOME/lib"

make
