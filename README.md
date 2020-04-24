# RTE+RRTMGP

This is the repository for RTE+RRTMGP, a set of codes for computing radiative fluxes in planetary atmospheres. RTE+RRTMGP is described in a [paper](https://doi.org/10.1029/2019MS001621) in [Journal of Advances in Modeling Earth Systems](http://james.agu.org).

RRTMGP uses a k-distribution to provide an optical description (absorption and possibly Rayleigh optical depth) of the gaseous atmosphere, along with the relevant source functions, on a pre-determined spectral grid given temperatures, pressures, and gas concentration. The k-distribution currently distributed with this package is applicable to the Earth's atmosphere under present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source functions. The fluxes are normally summarized or reduced via a user extensible class.

Example programs and documenation are evolving - please see examples/ in the repo and Wiki on the project's Github page. Suggestions are welcome. Meanwhile for questions please contact Robert Pincus and Eli Mlawer at rrtmgp@aer.com.

## Building the libraries.

This branch contains a version of the library for [ICON](https://code.mpimet.mpg.de/projects/iconpublic) with an [Autoconf](https://www.gnu.org/software/autoconf/)-based building system from [autoconf](https://github.com/RobertPincus/rte-rrtmgp/tree/autoconf).
The libraries can be configured, built, tested and installed following the standard workflow:

```shell
$ ./configure && make && make check && make install
```

The configure scripts supports the following use cases:
1. If you want to use a particular Fortran 2003 compiler with particular set of compiler flags:
    ```shell
    $ ./configure FC=<name of the compiler executable> FCFLAGS=<compiler flags>
    ```
2. If you want to use OpenACC kernels:
    ```shell
    $ ./configure --enable-openacc
    ```
3. If you want to build the [examples](#examples) by default (otherwise, they are built only for testing, i.e. when you call `make check`):
    ```shell
    $ ./configure --enable-examples
    ```
4. If you want to make sure that the tests are built and run (otherwise, the tests might be skipped):
    ```shell
    $ ./configure --enable-tests
    ```
5. If yout want to use [NetCDF Fortran library](https://www.unidata.ucar.edu/software/netcdf/docs-fortran/) (required for testing) from a non-standard directory:
    ```shell
    $ ./configure --with-netcdf-fortran=/path/to/netcdf-fortran
    ```

## Examples

Two examples are provided, one for clear skies and one including clouds. See the README file and codes in each directory for further information. 
