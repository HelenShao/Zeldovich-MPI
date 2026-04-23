# Standalone wisdom_toy: single-precision FFTW, same discovery + libs as parent Makefile default.
#
#   module load frameworks
#   module load fftw/3.3.10
#   make
#
# Convenience: make all-with-modules
#
# Override if needed (same variable roles as parent INCLUDES / FFTW_LIB_FLAG):
#   make FFTW_INCLUDES='-I/path' FFTW_LIB_FLAG='-L/path' FFTW_LIBS='$(FFTW_LIB_FLAG) -lfftw3f_omp -lfftw3f -lm'

CXX ?= g++
MPICXX ?= mpicxx
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra

# OpenMP: match parent (Intel -qopenmp vs GCC -fopenmp) so linking matches libfftw3f_omp.
ifeq ($(OPENMP_FLAGS),)
  OPENMP_FLAGS := $(shell echo "int main(){return 0;}" | $(CXX) -x c++ - -c -qopenmp -o /dev/null >/dev/null 2>&1 && echo -qopenmp || echo -fopenmp)
endif

# --- FFTW includes: same logic as ../Makefile (C_INCLUDE_PATH → /opt/aurora find → hardcoded) ---
ifneq ($(C_INCLUDE_PATH),)
  FFTW_INCLUDES := $(shell echo $(C_INCLUDE_PATH) | tr ':' '\n' | grep fftw | head -1 | sed 's|^|-I|')
  ifeq ($(strip $(FFTW_INCLUDES)),)
    FFTW_INCLUDES := $(shell find /opt/aurora -name "fftw3.h" 2>/dev/null | head -1 | xargs dirname | sed 's|^|-I|')
  endif
else
  FFTW_INCLUDES := -I/usr/local/fftw/gcc/openmpi-4.1.6/3.3.10/include
endif
ifeq ($(strip $(FFTW_INCLUDES)),)
  FFTW_INCLUDES := -I/usr/local/fftw/gcc/openmpi-4.1.6/3.3.10/include
endif

# --- FFTW lib dir: same logic as ../Makefile (LD_LIBRARY_PATH → /opt/aurora find → hardcoded) ---
ifneq ($(LD_LIBRARY_PATH),)
  FFTW_LIB_DIR := $(shell echo $(LD_LIBRARY_PATH) | tr ':' '\n' | grep fftw | head -1)
  ifeq ($(strip $(FFTW_LIB_DIR)),)
    FFTW_LIB_DIR := $(shell find /opt/aurora -name "libfftw3*.so" 2>/dev/null | head -1 | xargs dirname)
  endif
else
  FFTW_LIB_DIR :=
endif
ifeq ($(strip $(FFTW_LIB_DIR)),)
  FFTW_LIB_DIR := /usr/local/fftw/gcc/openmpi-4.1.6/3.3.10/lib
endif
FFTW_LIB_FLAG := -L$(FFTW_LIB_DIR)

# Single precision (default in parent ../Makefile): _omp + base float library
FFTW_LIBS := $(FFTW_LIB_FLAG) -lfftw3f_omp -lfftw3f -lm

ALL_CXXFLAGS := $(CXXFLAGS) $(OPENMP_FLAGS)

TARGET := wisdom_toy
SRC := wisdom_toy.cpp

SWEEP := fftw_2d_thread_sweep
SWEEP_SRC := fftw_2d_thread_sweep.cpp

SWEEP_MPI := fftw_2d_thread_sweep_mpi
SWEEP_MPI_SRC := fftw_2d_thread_sweep_mpi.cpp

.PHONY: all clean all-with-modules

all: $(TARGET) $(SWEEP) $(SWEEP_MPI)

all-with-modules:
	@chmod +x "$(CURDIR)/build_with_modules.sh" 2>/dev/null || true
	"$(CURDIR)/build_with_modules.sh" all

$(TARGET): $(SRC)
	$(CXX) $(ALL_CXXFLAGS) $(FFTW_INCLUDES) -o $@ $(SRC) $(FFTW_LIBS)

$(SWEEP): $(SWEEP_SRC)
	$(CXX) $(ALL_CXXFLAGS) $(FFTW_INCLUDES) -o $@ $(SWEEP_SRC) $(FFTW_LIBS)

$(SWEEP_MPI): $(SWEEP_MPI_SRC)
	$(MPICXX) $(ALL_CXXFLAGS) $(FFTW_INCLUDES) -o $@ $(SWEEP_MPI_SRC) $(FFTW_LIBS)

clean:
	rm -f $(TARGET) $(SWEEP) $(SWEEP_MPI)
