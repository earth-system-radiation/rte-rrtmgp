# RRTMGPXX testing
To run tests, execute the following:
```
cd ${HOME}
git clone https://github.com/mrnorman/YAKL.git
git clone https://github.com/E3SM-Project/rte-rrtmgp.git
cd rte-rrtmgp/cpp/test/build
source machine_environment_files/my_machine_file.sh
./get_data.sh
./cmakescript.sh 
make -j8 
make test
```

On some HPC machines, you'll need to clone rte-rrtmgp in a parallel filesystem directory instead, and you'll need to `make test` from within an interactive job.

