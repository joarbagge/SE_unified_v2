# The Spectral Ewald Unified package (version 2)

## Building

The code is built with CMake. To build, open a terminal in the
root directory and do

```
cd build
cmake ..
make
```

If CMake cannot find your Matlab installation, try

```
cmake .. -DMatlab_ROOT_DIR="/path/to/MATLAB/R2019a"
```

with the path replaced by the path to the Matlab installation.

## Demo scripts

Demo scripts are available in the ``SE{3,2,1,0}P`` directories.

## License

This package is released under the MIT License, see the [LICENSE](./LICENSE) file for details.

## Contributors

- Dag Lindbo
- Anna-Karin Tornberg
- Ludvig af Klinteberg
- Davood Saffar Shamshirgar
- Joar Bagge

## Research papers

A list of papers is found in the [PAPERS.md](./PAPERS.md) file.

## See also

See also the previous version of this package: <https://github.com/ludvigak/SE_unified>
