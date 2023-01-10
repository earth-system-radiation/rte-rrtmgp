### Contributing to RTE+RRTMGP

Thanks for considering making a contribution to RTE+RRTMGP.

The code in this repository is intended to work with compilers supporting the Fortran 2008 standard. It is also expected to run end-to-end on GPUs when compiled with OpenACC or OpenMP (though OpenMP is still unreliable). Commits are tested automatically against a range of compilers using Github Actions and also resources provided by the [Swiss Supercomputing Center](https://cscs.ch). The testing uses two general codes in `examples/`for which results are compared against existing implemetations,  and custom codes in `tests/` intended to excercise all code options.

##### Did you find a bug? 

Please file an issue on the [Github page](https://github.com/RobertPincus/rte-rrtmgp/issues). Please include a minimal reproducer of the bug it at all possible. 

##### Did you write a patch that fixes a bug?

Please fork this repository, branch from `develop`, make your changes, and open a Github [pull request](https://github.com/RobertPincus/rte-rrtmgp/pulls) against branch `develop`. 

##### Did you add functionality? 

Please fork this repository, branch from `develop`, make your changes, and open a Github [pull request](https://github.com/RobertPincus/rte-rrtmgp/pulls) against branch `develop`,  adding a new regression test or comparison against the reference in `tests/verification.py` or `tests/validation-plots.py` as appropriate. Add the test to the `tests` target in `tests/Makefile`.   

RTE+RRTMGP is intended to be a core that users can extend with custom code to suit their own needs. 