
# Location of the CUDA Toolkit
CUDA_PATH ?= /usr/local/cuda-8.0
TARGET_ARCH ?= x86_64
TARGET_SIZE := 64
HOST_COMPILER ?= g++
NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)

all:	run

clean:
	rm ProgressiveMesh

run:	main.cpp
	$(NVCC) main.cpp Render.cpp PM.cpp MeshLib_Core/Mesh.cpp MeshLib_Core/Vertex.cpp -o ProgressiveMesh -lGL -lGLU -lglut
