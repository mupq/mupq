# MuPQ build system

The build system consists of several platform-dependent (not part of this repo)
and platform independent modules. The build system assumes, that the `mupq`
repository is a submodule of the specific platform dependent repository, e.g.,
`pq{m4,m3,riscv}`. The build system is split into the following modules.

- `mupq/mk/config.mk`: The initial configuration script. It sets up all the
  basic platform independent compilation flags. It will also include the
  platform dependent configuration script named after the value of the
  `PLATFORM` variable. This configuration file also contains
  [mechanism](#variable-retention) for managing configuration variables that
  should trigger a `clean`.
- `mupq/mk/<PLATFORM>.mk`: This script sets up all the platform dependent
  compilation and linker flags. This file may contain a list of schemes to be
  skipped, meant for schemes that cannot run on the platform for any reason.
- `mupq/mk/rules.mk`: Contains all the pattern rules for compiling and linking
  source code.
- `mupq/mk/schemes.mk`: This script locates all the schemes, all their
  implementations and defines `make` targets all the objects, libraries and
  benchmark executables. The mechanism for locating schemes is defined
  [here](#locating-implementations).
- `mupq/mk/host-crypto.mk`: Defines a host-library for symmetric cryptographic
  algorithms that can be used by all the schemes. This host-library is platform
  independent and used to generate the test vectors.
- `mk/crypto.mk`: Defines the platform dependent symmetric crypto library. The
  library is compiled twice, once with- and once without profiling code enabled
  (used for profiling the time spent hashing / encrypting).

## User Configuration Variables

The build system can be parameterized by the user with the following variables:

- `PLATFORM=<yourplatform>` (required): The chosen target board/platform.
- `DEBUG=1`: Compile all code without optimization and with debug symbols.
- `OPT_SIZE=1`: Optimize all code for size (otherwise the default is `-O3`).
- `LTO=1`: Enable link-time optimization.
- `AIO=1`: Use all-in-one compilation of schemes, i.e. pass all sources instead
  of compiled modules to the linking step (this can, in some cases, be faster
  than link-time optimization).

Each platform may also define further variables.

## Output folders

The resulting binaries are placed in the following folders:

- `obj`: The compiled objects and libraries.
- `obj-host`: The compiled host-objects and host-libraries.
- `elf`: The linked `elf` executables of all the benchmarks.
- `bin`: A flashable binary file (if required by the target platform).
- `bin-host`: The linked executables for the host, e.g., test vector generators.

## Variable Retention

The build system will retain the value of certain variables by generating a
`obj/.config.mk` that contains the values of all retained variables. The
variables that will be retained are defined by the `RETAINED_VARS` variable. An
error will be triggered promting the user to run `make clean`, if any of these
variables appear to be changed from the initial configuration value. This is
meant for configuration variables that affect the compilation result, such as
optimization flags.
  
## Locating Implementations

Each time `make` is run, a set of folders will be searched for implementations
to define targets for libraries and benchmark executables. Schemes are meant to
located in the following folders:

- `crypto_{kem,sign}/`: Platform dependent KEM/signature schemes.
- `mupq/crypto_{kem,sign}/`: Platform independent KEM/signature schemes.
- `mupq/pqclean/crypto_{kem,sign}/`: Platform independent KEM/signature schemes
  from the PQClean project.

Each scheme is placed in a subfolder named after a scheme, with each
implementation in its own folder, i.e.,
`crypto_{kem,sign}/<scheme>/<implementation>/<sourcefiles.{c,h,s}>`.

The `mupq/mk/schemes.mk` will automatically find these folders and define
targets for each implementation. The name of each implementation is simply
derived from its path, with `/` replaced by `_`. The following targets are
defined:

- A library `obj/lib<schemename>.a`
- A self test `elf/<schemename>_test.elf`
- A speed test `elf/<schemename>_speed.elf`
- A hashing profiling test `elf/<schemename>_hashing.elf`
- A stack size test `elf/<schemename>_stack.elf`
- A test vector test `elf/<schemename>_testvector.elf`
- If the implementation is platform independent, a test vector generator
  `bin-host/<schemename>_testvectors` will be defined.
- `bin/*.bin` / `bin/*.hex` files are generated automatically for each `elf`
  file, if required by the flashing tools for the platform.

The locating mechanism of the script utilizes the capability of `make` to
[(re-)make included makefiles](https://www.gnu.org/software/make/manual/html_node/Remaking-Makefiles.html).
The script will automatically generate a `obj/.schemes.mk` files, that contains
the two variables `KEM_SCHEMES` and `SIGN_SCHEMES`, each of which is just a list
of paths to all implementations. If this file does not exist, it will be
generated by two calls to the `find` utility that will enumerate all
implementations within the search paths mentioned above. The `obj/.schemes.mk`
is furthermore marked "phony", triggering `make` to re-generate it everytime it
is called.

This location mechanism will be skipped, if the `IMPLEMENTATION_PATH` variable
is defined (all `make` invocations by the python benchmarking scripts will
automatically include this variable). In this case, only the single
implementation pointed to by this variable will be compiled.
