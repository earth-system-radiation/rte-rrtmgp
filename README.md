# RTE+RRTMGP

This is the repository for RTE+RRTMGP, a set of codes for computing radiative fluxes in planetary atmospheres. RTE+RRTMGP is described in a [paper](https://doi.org/10.1029/2019MS001621) in [Journal of Advances in Modeling Earth Systems](http://james.agu.org).

RRTMGP uses a k-distribution to provide an optical description (absorption and possibly Rayleigh optical depth) of the gaseous atmosphere, along with the relevant source functions, on a pre-determined spectral grid given temperatures, pressures, and gas concentration. The k-distribution currently distributed with this package is applicable to the Earth's atmosphere under present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source functions. The fluxes are normally summarized or reduced via a user extensible class.

Example programs and documentation are evolving - please see examples/ in the repo and Wiki on the project's Github page. Suggestions are welcome. Meanwhile for questions please contact Robert Pincus and Eli Mlawer at rrtmgp@aer.com.

## Recent changes

1. The default method for solution for longwave problems that include scattering has been changed from 2-stream methods to a re-scaled and refined no-scattering calculation following [Tang et al. 2018](https://doi.org/10.1175/JAS-D-18-0014.1).
2. In RRTMGP gas optics, the spectrally-resolved solar source function in can be adjusted by specifying the total solar irradiance (`gas_optics%set_tsi(tsi)`) and/or the facular and sunspot indicies (`gas_optics%set_solar_variability(mg_index, sb_index, tsi)`)from the [NRLSSI2 model of solar variability](http://doi.org/10.1175/BAMS-D-14-00265.1).  
3. `rte_lw()` now includes optional arguments for computing the Jacobian (derivative) of broadband flux with respect to changes in surface temperature. In calculations neglecting scattering only the Jacobian of upwelling flux is computed. When using re-scaling to account for scattering the Jacobians of both up- and downwelling flux are computed.

Relative to commit `69d36c9` to `master` on Apr 20, 2020, the required arguments to both the longwave and shortwave versions of `ty_gas_optics_rrtmgp%load()`have changed.


## Building the libraries.

This branch contains an [Autoconf](https://www.gnu.org/software/autoconf/)-based building system for [master](https://github.com/RobertPincus/rte-rrtmgp/tree/master).
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

Two [examples](./examples) are provided, one for [clear skies](./examples/rfmip-clear-sky) and one [including clouds](./examples/all-sky). See the README file and codes in each directory for further information.
