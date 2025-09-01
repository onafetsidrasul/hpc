# Compilers
CC = gcc
MPICC = mpicc

# Build options
OPENMP_SCHEDULE ?= static      # [static|dynamic|guided|auto|0-3]
BUILD_TYPE ?= release          # [release|debug|profile]
ARCH ?= native                 # CPU architecture (native for auto-detection)
OMP_CHUNK_SIZE=1000

# Map OpenMP schedule types to integer flags
ifeq ($(OPENMP_SCHEDULE), static) 
  SCHEDULE_FLAG = 0
  OMP_CHUNK_SIZE = 1
else ifeq ($(OPENMP_SCHEDULE), dynamic)
  SCHEDULE_FLAG = 1
else ifeq ($(OPENMP_SCHEDULE), guided)
  SCHEDULE_FLAG = 2
else ifeq ($(OPENMP_SCHEDULE), auto)
  SCHEDULE_FLAG = 3
else
  SCHEDULE_FLAG = $(OPENMP_SCHEDULE)
endif

# Base compiler flags
BASE_FLAGS = -Wall -Wextra -I./include -fopenmp -march=$(ARCH)
MPI_BASE_FLAGS = $(BASE_FLAGS)

# Build type specific flags
ifeq ($(BUILD_TYPE), debug)
  CFLAGS = $(BASE_FLAGS) -g -O0 -DDEBUG -fno-omit-frame-pointer
  MPI_FLAGS = $(MPI_BASE_FLAGS) -g -O0 -DDEBUG -fno-omit-frame-pointer
else ifeq ($(BUILD_TYPE), profile)
  CFLAGS = $(BASE_FLAGS) -g -O3 -pg -DPROFILE
  MPI_FLAGS = $(MPI_BASE_FLAGS) -g -O3 -pg -DPROFILE
else
  # Release flags
  CFLAGS = $(BASE_FLAGS) -O3 -DNDEBUG -flto
  MPI_FLAGS = $(MPI_BASE_FLAGS) -O3 -DNDEBUG -flto
endif

# Add OpenMP
CFLAGS += -DOPENMP_SCHEDULE=$(SCHEDULE_FLAG) -DOMP_CHUNK_SIZE=$(OMP_CHUNK_SIZE)
MPI_FLAGS += -DOPENMP_SCHEDULE=$(SCHEDULE_FLAG) 


# Non-OpenMP builds
NO_OMP_FLAGS = $(filter-out -fopenmp,$(BASE_FLAGS)) -I./include

# Source directory
SRC_DIR = src

# Source files
SRC_SERIAL       = $(SRC_DIR)/stencil_template_serial.c
SRC_SERIAL_FR    = $(SRC_DIR)/stencil_template_serial_original.c
SRC_PARALLEL     = $(SRC_DIR)/stencil_template_parallel.c

# Object files
OBJ_SERIAL       = $(SRC_SERIAL:.c=.o)
OBJ_SERIAL_FR    = $(SRC_SERIAL_FR:.c=.o)
OBJ_PARALLEL     = $(SRC_PARALLEL:.c=.o)

# Executables
EXEC_SERIAL      = stencil_serial
EXEC_SERIAL_FR   = stencil_serial_fr
EXEC_PARALLEL    = stencil_parallel

# Default target: build everything in release mode
all: $(EXEC_SERIAL) $(EXEC_SERIAL_FR) $(EXEC_PARALLEL)

# Build only serial versions
openmp: $(EXEC_SERIAL) $(EXEC_SERIAL_FR)

# Serial with OpenMP
$(EXEC_SERIAL): $(OBJ_SERIAL)
	$(CC) $(CFLAGS) -o $@ $^

# Fallback serial without OpenMP
$(EXEC_SERIAL_FR): $(OBJ_SERIAL_FR)
	$(CC) $(NO_OMP_FLAGS) -o $@ $^

# Parallel (MPI + OpenMP)
$(EXEC_PARALLEL): $(OBJ_PARALLEL)
	$(MPICC) $(MPI_FLAGS) -o $@ $^

# Individual object compilation
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean targets
clean:
	rm -f $(SRC_DIR)/*.o $(EXEC_SERIAL) $(EXEC_SERIAL_FR) $(EXEC_PARALLEL) 

allclean:
	rm -f $(SRC_DIR)/*.o $(EXEC_SERIAL) $(EXEC_SERIAL_FR) $(EXEC_PARALLEL) *.txt *.bin *.csv *.out gmon.out

# Help target
help:
	@echo "Build options:"
	@echo "  make BUILD_TYPE=debug    - Debug build (-g -O0)"
	@echo "  make BUILD_TYPE=profile  - Profiling build (-pg)"
	@echo "  make BUILD_TYPE=release  - Optimized release (default)"
	@echo ""
	@echo "OpenMP scheduling:"
	@echo "  make OPENMP_SCHEDULE=static|dynamic|guided|auto"
	@echo ""
	@echo "Architecture:"
	@echo "  make ARCH=native         - Optimize for current CPU (default)"
	@echo "  make ARCH=x86-64-v3      - Specific architecture"

.PHONY: all clean allclean openmp help
