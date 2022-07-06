include etc/cc4x.mk

.DEFAULT_GOAL := install
.PHONY: install clean

SRC_PATH          = ./src
BUILD_PATH        = build
OBJ_PATH          = $(BUILD_PATH)/obj
BIN_PATH          = $(BUILD_PATH)/bin
BIN_PROGRAMS      = $(BIN_PATH)/${EXE}

SRC := $(wildcard $(SRC_PATH)/*.cxx)
OBJ_FILES = $(patsubst $(SRC_PATH)/%.cxx,$(OBJ_PATH)/%.o,$(SRC))

clean:
	rm -rf $(BUILD_PATH)

EXE := cc4x
install: $(BIN_PROGRAMS)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.cxx
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_PATH)/${EXE}: ${OBJ_FILES}
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) ${OBJ_FILES} ${LIBS} -o $@


