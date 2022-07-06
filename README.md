# Scope of the project 

`cc4x` is a canoncial coupled cluster calculations code. It uses tensor
contraction libraries for the occurrent tensor operations. Currently, only CTF
(cyclops) is supported.  However, adding other contraction libraries is
straightforward (see TensorBackend.md). The offered functionalities are
limited! Please check out [cc4s](https://github.com/cc4s/cc4s) for general
survey. 


# Software dependencies

You will need [yaml-cpp](https://github.com/jbeder/yaml-cpp/) and
[cyclops](https://github.com/cyclops-community/ctf). Both can be found on
github.

# Data dependencies

When running the uniform electron gas, a well known model system in solid state
physics, no input files are required. Reasonal input parameters can be found in
the work of Morales and co-workers
[afqmc-ueg](https://aip.scitation.org/doi/pdf/10.1063/1.5109572).  Furthermore,
the input files from `cc4s` can be used (with tiny modifications in the yaml).

# How to install

Go to the file etc/cc4x.mk and define the ```CTF_PATH``` and ```YAML_PATH```.
Then go to root directory and type ```make```. This build system assumes gcc
compiler with openmpi and OpenBlas.


# How to get started

You might want to start running the BN test provided in the ```test```
directory.  You should be able to reproduce the CCSD energy, calculated with
`cc4s`.  Simply run `cc4x`: ```$cc4x -i 40``` If you want to run some (large
scale) calculations on the uniform electron gas.  For instance ```$cc4x  -o 27
-v 224 -r 1``` Note that the value of ```r``` will have no effect on the
performance of the calculation (simply choose a positiv value; 1.0 is fine).
Have a look at
[cc4s-manual](https://manuals.cc4s.org/user-manual/performance/performance.html)
for some large scale calculations of `cc4s`.
