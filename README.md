# Point Set Analyser

## Author
Jean-Yves Franceschi <jean-yves.franceschi@ens-lyon.org>

## About

## Requirements
The following dependencies are required:
 - a C++ compiler supporting C++11 standard;
 - cmake;
 - Qt 5.6 or above;
 - CGAL 4.7 or above.
In order to speed up the computations, the C++ compiler should provide OpenMP's support.

## Compilation
To compile you need to (on Linux systems):
 - in the root folder of the project, create a ```build``` folder;
 - in the build folder, run ```cmake ..``` then ```make```.

## Usage
The executable is stored in the ```bin``` folder of the root directory of the project. It just needs to be executed without any options or arguments.

The loaded point sets must be text files of the following form: on a first line it shows the number of points N, and on the following N lines it displays a pair of doubles between 0 and 1 representing the coordinated of a point of the set.

For the moment only PGM images can be loaded in the viewer.

## Doxygen documentation
If Doxygen is available, then ```make doc``` in ```build``` produces a documentation to be found in ```doc```.
