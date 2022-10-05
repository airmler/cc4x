
CXX=CC
CXXFLAGS := -std=c++11 -O3 -fopenmp



CTF_PATH = /users/irmleran/src/cc4x/ctf
YAML_PATH = /users/irmleran/src/cc4s/extern/build/lumi-gcc/yaml-cpp/c9460110e072df84b7dee3eb651f2ec5df75fb18


# BACKENDS
CXXFLAGS += -DCTF_BACKEND -DCONTRACTION_COLLECTOR


CXXFLAGS += -I${CTF_PATH}/include
CXXFLAGS += -I. -Isrc/
CXXFLAGS += -I${YAML_PATH}/include

LIBS += -L${CTF_PATH}/lib -lctf
#LIBS += -L/opt/OpenBLAS/lib/ -lopenblas -fopenmp -lscalapack -lgfortran
LIBS += -lyaml-cpp -L${YAML_PATH}/lib64

