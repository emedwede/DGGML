# Update: Please Read
The build instructions are out of data and will be updated soon! Thanks for your patience. 

# DGGML: The Dynamical Graph Grammar Modeling Library

DGGML is a modeling library for developing simulations of spatially embedded dynamical graph grammars. It's a follow up to [CajeteCMA](https://github.com/emedwede/CajeteCMA) - a Graph Grammar Simulator for the Cortical Microtubule Array (CMA), which has a serial implementation for a graph grammar simulator featured in the paper titled, ["Approximate Simulation of Cortical Microtubule Models using Dynamical Graph Grammars"](https://dx.doi.org/10.1088/1478-3975/acdbfb). The graph library used is YAGL: Yet another Graph Library.


# Requirements and Build Instructions

## Linux
The minimum CXX version required to compile is C++17. The code has been tested on Pop!_OS 20.04 LTS, but most versions of linux should be able to compile and run the code. We also use the state of the art ODE solver [SUNDIALS](https://github.com/LLNL/sundials), and a build of the latest version is required to be provided in the build script for CMAKE to find. Every other requirement is packaged in. The build instructions are fairly simple, but require a little bit of work. First build [SUNDIALS](https://github.com/LLNL/sundials), then read and modify the builder script in the scripts directory and simply run ./scripts/builder.sh The build will be placed out of source wherever the user defines.

## MacOS
To get the code to run on the new ARM based Apple silicon(poorly named M1/M2/pro/max/etc), a few commands in homebrew will get you going. It'll probably work on older intel based chips too. If you haven't already, check out [homebrew](https://brew.sh/) Why not live a little? Once you have brew installed, follow these instructions:

1) brew install sundials
2) brew install brew install gcc@12
3) brew install cmake
4) clone this repository to wherever you want
6) cd into the repo, type `mkdir build && cd build`
7) from the build directory set the compilers and configure the build with `CC=gcc-12 CXX=g++12 cmake ..`
8) Then to build `make -j4`

Please not that the clang compiler is currently untested, so use at your discretion.
