BIN = growth
CXX = clang++
NVCC = nvcc
CUDA_DIR = /usr/local/cuda
NVCC_CLANG = clang++

OBJ =	build/growth/main.o \
		build/growth/Knots.o \
		build/growth/RD.o \
		build/growth/Flow.o \
		build/growth/growth.o \
		build/growth/WriteData.o \
		build/Assembler.o \
		build/LevelSet.o \
		build/TetraMesh.o

CUDA_OBJ = build/EuclideanDistance.o

CFLAGS = -DCFISH_GROWTH -D_FORTIFY_SOURCE=2 -DNDEBUG
CFLAGS += -DCGAL_LINKED_WITH_TBB -DCGAL_USE_GMP -DCGAL_USE_MPFR
CFLAGS += -std=c++14 -fopenmp -march=native -O3 -Wall -Wno-deprecated-declarations
CFLAGS += -I./src -I./src/$(BIN) -I$(CUDA_DIR)/include

# OS-specific flags; currently only macOS and Linux (default) recognized.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
CFLAGS += -I/opt/local/include/
else
CFLAGS +=
endif

NVCC_FLAGS = -std=c++14 -O3 -Xcompiler "-Wall -march=native -fopenmp" \
		     -ccbin $(NVCC_CLANG) -Xptxas -O3,-v

LDFLAGS = -lgmp -fopenmp -ltbb -ltbbmalloc -lCGAL
LDFLAGS += -lcudart -Wl,-rpath,$(CUDA_DIR)/lib
ifeq ($(UNAME_S),Darwin)
LDFLAGS += -L$(CUDA_DIR)/lib -L/opt/local/lib
else
LDFLAGS += -L$(CUDA_DIR)/targets/x86_64-linux/lib
endif

all: bin/$(BIN)

build/$(BIN):
	@mkdir -p build/$(BIN)

build/$(BIN)/%.o : src/%.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

build/$(BIN)/%.o : src/$(BIN)/%.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

build/%.o : src/fem/%.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

build/%.o : src/common/%.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

build/%.o : src/fem/%.cu
	$(NVCC) $(NVCC_FLAGS) -c -o $@ $<

bin/$(BIN): build/$(BIN) $(OBJ) $(CUDA_OBJ)
	$(CXX) -o bin/$(BIN) $(OBJ) $(CUDA_OBJ) $(LDFLAGS)

.PHONY: clean

clean:
	rm -rf build/* bin/*

install:
	cp $(BIN) ../bin/

