# Scope of the project 

`cc4x` is a canoncial coupled cluster calculations code. It uses tensor
contraction libraries for the occurrent tensor operations. Currently, only CTF
(cyclops) is supported.  However, adding other contraction libraries is
straightforward (see TensorBackend.md). The offered functionalities are
limited! Please check out [cc4s](https://github.com/cc4s/cc4s) for general
survey. 


# Software dependencies

You will need [yaml-cpp](https://github.com/jbeder/yaml-cpp/) and
[cyclops](https://github.com/cyclops-community/ctf) which can be found on
github.

# Data dependencies

When running the uniform electron gas, a well known model system in solid state
physics, no input files are required. Reasonal input parameters can be found in
the work of Morales and co-workers
[afqmc-ueg](https://aip.scitation.org/doi/pdf/10.1063/1.5109572).  Furthermore,
the input files from `cc4s` can be used (with tiny modifications in the yaml).

# How to install

Go to the file etc/cc4x.mk and define the ```CTF_PATH``` and ```YAML_PATH```.
Then go to ```src``` and type ```make```. This build system assumes gcc
compiler with openmpi and OpenBlas.
