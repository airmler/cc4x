
CXX=mpic++
#CXXFLAGS := -std=c++11 -O3 -fmax-errors=1
CXXFLAGS := -std=c++11 -O0 -g -fmax-errors=1
#CXXFLAGS := -std=c++11 -O0 -g -fmax-errors=1 -Wall -pedantic



CTF_PATH = /home/irmler/Programs/cc4x/ctf
YAML_PATH = /home/irmler/Programs/cc4s_release/extern/build/gcc-oblas-ompi/yaml-cpp/c9460110e072df84b7dee3eb651f2ec5df75fb18


# BACKENDS
CXXFLAGS += -DIBST_BACKEND
#CXXFLAGS += -DCTF_BACKEND


CXXFLAGS += -I${CTF_PATH}/include
CXXFLAGS += -Iibst/src
CXXFLAGS += -I. -Isrc/
CXXFLAGS += -I${YAML_PATH}/include

LIBS += -L${CTF_PATH}/lib -lctf
LIBS += -L/opt/OpenBLAS/lib/ -lopenblas -fopenmp -lscalapack -lgfortran
LIBS += -lyaml-cpp -L${YAML_PATH}/lib

