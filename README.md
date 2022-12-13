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
