### Contributing to RTE+RRTMGP

Thanks for considering making a contribution to RTE+RRTMGP.

The code in this repository is intended to work with compilers supporting the Fortran 2008 standard. It is also expected to run end-to-end on GPUs when compiled with OpenACC. Commits are tested with [Travis](https://travis-ci.com) against gfortran versions 8 and 9 and against various versions > 19.9 of the PGI compiler using resources provided by the [Swiss Supercomputing Center](https://cscs.ch). The testing uses two general codes in `examples/`for which results are compared against existing implemetations,  and custom codes in tests/ intended to excercise all code options.

##### Did you find a bug? 

Please file an issue on the [Github page](https://github.com/RobertPincus/rte-rrtmgp/issues). 

##### Did you write a patch that fixes a bug?

Please fork this repository, branch from `develop`, make your changes, and open a Github [pull request](https://github.com/RobertPincus/rte-rrtmgp/pulls) against branch `develop`. 

##### Did you add functionality? 

Please fork this repository, branch from `develop`, make your changes, and open a Github [pull request](https://github.com/RobertPincus/rte-rrtmgp/pulls) against branch `develop`,  adding a new regression test or comparison against the reference in `tests/verification.py` or `tests/validation-plots.py` as appropriate.  

RTE+RRTMGP is intended to be a core that users can extend with custom code to suit their own needs. 