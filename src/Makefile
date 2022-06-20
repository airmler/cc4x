.DEFAULT_GOAL := install
.PHONY: install clean

BUILD_PATH        = build
OBJ_PATH          = $(BUILD_PATH)/obj
BIN_PATH          = $(BUILD_PATH)/bin

SRC := $(wildcard *.cxx)
OBJ_FILES = $(patsubst %.cxx,%.o,$(SRC))

clean:
	rm *.o
	rm cc4x

CXX=mpic++
CXXFLAGS := -std=c++11 -g -O0


# BACKENDS
CXXFLAGS += -DCTFBS_BACKEND


CTF_PATH = /home/irmler/uni/kccsd/cc4x/ctf
YAML_PATH = /home/irmler/uni/kccsd/cc4x/yaml-cpp/build/install
CXXFLAGS += -I${CTF_PATH}/include
CXXFLAGS += -I.
CXXFLAGS += -I${YAML_PATH}/include

LIBS += -L${CTF_PATH}/lib -lctf
LIBS += -lopenblas -fopenmp
LIBS += -lyaml-cpp -L${YAML_PATH}/lib
EXE := cc4x

install: ${EXE}

%.o: %.cxx
	$(CXX) -c $(CXXFLAGS) $<

${EXE}: ${OBJ_FILES}
	$(CXX) $(CXXFLAGS) ${OBJ_FILES} ${LIBS} -o $@


