# TPMSproject
This project aims at implementing a framework to create different types of Triply Periodic Minimal Surfaces (TPMS).

## Features
This code is written in [C++](https://isocpp.org/) and uses the [VTK](https://vtk.org/) tool to create the TPMS geometry and decimate the output surface mesh.
Both *sheet* and *solid* TPMS are implemented, for several topologies.

The code has been written based on another GITHUB project [TPMSGenerator](https://github.com/EGajetti/TPMSgenerator).

## TPMS topologies
Several TPMS types, or topologies, are implemented: 
- Primitive
- Gyroid
- Diamond
- IWP
- Fischer-Koch-S
- F-RD

## Install
In order to use this code, VTK should be installed, as reported [here](https://gitlab.kitware.com/vtk/vtk/-/blob/v9.2.6/Documentation/dev/build.md#building-vtk), possibly the latest version (VTK9). If you are using a previous VTK version, the CCMakeLists might be slightly different, as reported [here](https://vtk.org/doc/nightly/html/md__builds_gitlab-kitware-sciviz-ci_Documentation_Doxygen_ModuleMigration.html).

When having VTK installed on your pc, clone the repository in your preferred folder:

$ git clone https://github.com/MAHTEP/TPMSproject.git

Then create a new folder *build* and enter it:

$ mkdir build
$ cd build

Now you need to compile the project, using:

$ cmake .. -GNinja
$ cmake --build .
