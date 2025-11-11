# sharp2DQR

This is an implementation of sharp2DQR for counting the number of Skolem functions of a given 2-DQBF.

The theoretical study of #2-DQBF that forms the basis of sharp2DQR is published in the paper:

  Model Counting for Dependency Quantified Boolean Formulas\
  Long-Hin Fung, Che Cheng, Jie-Hong Roland Jiang, Friedrich Slivovsky, Tony Tan\
  To appear in the Proceedings of AAAI 2026

## Building this project
- Install [`vcpkg`](https://github.com/microsoft/vcpkg)
- Install `cxxopts` with `{VCPKG_ROOT}/vcpkg install cxxopts`
- Run the following commands

```
mkdir build
cd build
cmake -DCMAKE_TOOLCHAIN_FILE={VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake ..
make
```
- (MacOS) Please install [`homebrew`](https://brew.sh) and run `brew install coreutils`

## Usage

```sharp2DQR
Usage:
  sharp2DQR [OPTION...]

  -i, --input arg   Input File
  -v, --verbose     Verbose output
  -h, --help        Print usage
```
