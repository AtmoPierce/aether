SHELL := /bin/bash

CUDA_MODULE ?= cuda/12.8.1
MODULE_INIT ?= /apps/profiles/modules_asax.sh.dyn
CUDA_HOME ?=
CUDA_GENCODE ?=

BUILD_DIR ?= build
OBJ_DIR := $(BUILD_DIR)/obj
PTX_DIR := $(BUILD_DIR)/ptx
LIB_NAME ?= libfpx_cuda_kernels.a
LIB_TARGET := $(BUILD_DIR)/$(LIB_NAME)
PERF_SRC := cuda_perf_main.cu
PERF_TARGET := $(BUILD_DIR)/fpx_cuda_perf
PRECISION_SRC := cuda_precision_main.cu
PRECISION_TARGET := $(BUILD_DIR)/fpx_cuda_precision
MAIN_SOURCES := $(PERF_SRC) $(PRECISION_SRC)

CU_SOURCES := $(sort $(filter-out $(MAIN_SOURCES),$(wildcard *.cu)))
OBJECTS := $(patsubst %.cu,$(OBJ_DIR)/%.o,$(CU_SOURCES))
PTX_FILES := $(patsubst %.cu,$(PTX_DIR)/%.ptx,$(CU_SOURCES))
PERF_OBJ := $(OBJ_DIR)/$(PERF_SRC:.cu=.o)
PRECISION_OBJ := $(OBJ_DIR)/$(PRECISION_SRC:.cu=.o)

NVCC_STD ?= -std=c++17
OPTFLAGS ?= -O3
PICFLAGS ?= -Xcompiler -fPIC
WARNFLAGS ?= -Xcudafe --display_error_number
DEFINES ?= -DFPX_CUDA_USE_FMA=1
INCLUDES ?=
NVCCFLAGS := $(strip $(NVCC_STD) $(OPTFLAGS) $(PICFLAGS) $(WARNFLAGS) $(DEFINES) $(CUDA_GENCODE) $(INCLUDES))

AR ?= ar
ARFLAGS ?= rcs

ifeq ($(strip $(CUDA_HOME)),)
NVCC = source $(MODULE_INIT) >/dev/null 2>&1 && module load $(CUDA_MODULE) >/dev/null 2>&1 && nvcc
else
NVCC = $(CUDA_HOME)/bin/nvcc
endif

.PHONY: all bench precision objects ptx clean print-config help

all: $(LIB_TARGET)

bench: $(PERF_TARGET)

precision: $(PRECISION_TARGET)

objects: $(OBJECTS)

ptx: $(PTX_FILES)

$(BUILD_DIR) $(OBJ_DIR) $(PTX_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: %.cu | $(OBJ_DIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(PTX_DIR)/%.ptx: %.cu | $(PTX_DIR)
	$(NVCC) $(NVCCFLAGS) -ptx $< -o $@

$(LIB_TARGET): $(BUILD_DIR) $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

$(PERF_TARGET): $(BUILD_DIR) $(OBJECTS) $(PERF_OBJ)
	$(NVCC) $(NVCCFLAGS) $(PERF_OBJ) $(OBJECTS) -o $@

$(PRECISION_TARGET): $(BUILD_DIR) $(OBJECTS) $(PRECISION_OBJ)
	$(NVCC) $(NVCCFLAGS) $(PRECISION_OBJ) $(OBJECTS) -o $@

clean:
	rm -rf $(BUILD_DIR)

print-config:
	@printf 'CUDA_MODULE=%s\n' "$(CUDA_MODULE)"
	@printf 'MODULE_INIT=%s\n' "$(MODULE_INIT)"
	@printf 'CUDA_HOME=%s\n' "$(CUDA_HOME)"
	@printf 'CUDA_GENCODE=%s\n' "$(if $(strip $(CUDA_GENCODE)),$(CUDA_GENCODE),<nvcc default>)"
	@printf 'NVCC=%s\n' "$(NVCC)"
	@printf 'BUILD_DIR=%s\n' "$(BUILD_DIR)"
	@printf 'PERF_TARGET=%s\n' "$(PERF_TARGET)"
	@printf 'PRECISION_TARGET=%s\n' "$(PRECISION_TARGET)"
	@printf 'SOURCES=%s\n' "$(CU_SOURCES)"

help:
	@printf '%s\n' \
		'Targets:' \
		'  make / make all      Build every .cu file and archive them into $(LIB_TARGET)' \
		'  make bench           Build the standalone CUDA performance benchmark $(PERF_TARGET)' \
		'  make precision       Build the standalone CUDA precision benchmark $(PRECISION_TARGET)' \
		'  make objects         Build only object files in $(OBJ_DIR)' \
		'  make ptx             Emit PTX files into $(PTX_DIR)' \
		'  make clean           Remove build outputs' \
		'' \
		'Overrides:' \
		'  make CUDA_MODULE=cuda/12.8.1' \
		'  make CUDA_HOME=/usr/local/cuda-12.8' \
		'  make CUDA_GENCODE="-gencode arch=compute_90a,code=sm_90a"'