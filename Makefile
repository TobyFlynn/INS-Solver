MAKEFILES_DIR != dirname $(realpath \
  $(word $(words $(MAKEFILE_LIST)), $(MAKEFILE_LIST)))
ROOT_DIR != realpath $(MAKEFILES_DIR)

INC := -I$(ROOT_DIR)/include
BIN := $(ROOT_DIR)/bin
OBJ := $(ROOT_DIR)/obj
2D_SEQ_OBJ_DIR := $(OBJ)/2d_seq
2D_OMP_OBJ_DIR := $(OBJ)/2d_omp
2D_CUDA_OBJ_DIR := $(OBJ)/2d_cuda
2D_HIP_OBJ_DIR := $(OBJ)/2d_hip
2D_MPI_SEQ_OBJ_DIR := $(OBJ)/2d_mpi_seq
2D_MPI_OMP_OBJ_DIR := $(OBJ)/2d_mpi_omp
2D_MPI_CUDA_OBJ_DIR := $(OBJ)/2d_mpi_cuda
2D_MPI_HIP_OBJ_DIR := $(OBJ)/2d_mpi_hip
3D_SEQ_OBJ_DIR := $(OBJ)/3d_seq
3D_OMP_OBJ_DIR := $(OBJ)/3d_omp
3D_CUDA_OBJ_DIR := $(OBJ)/3d_cuda
3D_HIP_OBJ_DIR := $(OBJ)/3d_hip
3D_MPI_SEQ_OBJ_DIR := $(OBJ)/3d_mpi_seq
3D_MPI_OMP_OBJ_DIR := $(OBJ)/3d_mpi_omp
3D_MPI_CUDA_OBJ_DIR := $(OBJ)/3d_mpi_cuda
3D_MPI_HIP_OBJ_DIR := $(OBJ)/3d_mpi_hip

include config.mk

OP2_INC = -I$(OP2_DIR)/include
OPENBLAS_INC = -I$(OPENBLAS_DIR)/include
ARMA_INC = -I$(ARMA_DIR)/include
PETSC_INC = -I$(PETSC_DIR)/include
INIPP_INC = -I$(INIPP_DIR)
MPI_INC = -I$(MPI_DIR)/include
OP2_DG_TOOLKIT_INC = -I$(OP2_DG_TOOLKIT_DIR)/include

OPENBLAS_LIB = -L$(OPENBLAS_DIR)/lib -lopenblas
PETSC_LIB = -L$(PETSC_DIR)/lib -lpetsc
ARMA_LIB = -L$(ARMA_DIR)/lib64 -larmadillo
HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5
COMMON_LIBS = $(OPENBLAS_LIB) $(PETSC_LIB) $(ARMA_LIB) $(HDF5_LIB)

PARTITION_LIB = -L$(PARMETIS_DIR)/lib -lparmetis -lmetis -lGKlib

OP2_OPENMP_LIBS = -L$(OP2_DIR)/lib -lop2_openmp -lop2_hdf5
OP2_MPI_LIBS = -L$(OP2_DIR)/lib -lop2_mpi

# Common compile flags
BASE_FLAGS := -g -O3
SEQ_CPU_COMPILER_FLAGS := $(BASE_FLAGS)
OMP_CPU_COMPILER_FLAGS := $(BASE_FLAGS) $(OPENMP_FLAG)
CUDA_COMPILER_FLAGS := $(BASE_FLAGS) $(OPENMP_FLAG)
NVCC_FLAGS := $(BASE_FLAGS) -rdc=true -gencode arch=compute_$(CUDA_ARCH),code=sm_$(CUDA_ARCH) -Xcompiler $(OPENMP_FLAG) $(MPI_INC)
HIP_FLAGS := $(BASE_FLAGS) -fgpu-rdc --offload-arch=$(HIP_ARCH) $(OPENMP_FLAG) $(MPI_INC)
COMMON_COMPILE_DEFS := -DDG_ORDER=$(ORDER) -DDG_COL_MAJ -DMAX_CONST_SIZE=1024 -DOP2_PARTITIONER=$(PART_LIB_NAME)
COMMON_COMPILE_DEFS_2D := $(COMMON_COMPILE_DEFS) -DDG_DIM=2
COMMON_COMPILE_DEFS_3D := $(COMMON_COMPILE_DEFS) -DDG_DIM=3
HIP_COMPILE_DEFS := -DOP2_DG_CUDA -DDG_OP2_SOA
CUDA_COMPILE_DEFS := -DOP2_DG_CUDA -DDG_OP2_SOA
MPI_COMPILER_DEFS := -DDG_MPI -DINS_MPI
INS_INC := $(ARMA_INC) $(INC) $(OP2_INC) $(OPENBLAS_INC) $(PETSC_INC) $(OP2_DG_TOOLKIT_INC) $(INIPP_INC)

all: cpu

cpu: $(BIN)/ins2d_openmp $(BIN)/ins2d_mpi $(BIN)/ins2d_mpi_openmp \
	$(BIN)/ins3d_openmp $(BIN)/ins3d_mpi $(BIN)/ins3d_mpi_openmp

clean:
	-rm -rf $(OBJ)
	-rm -rf $(BIN)

$(BIN):
	@mkdir -p $@

$(OBJ):
	@mkdir -p $@

$(2D_SEQ_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(2D_OMP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(2D_CUDA_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(2D_HIP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(2D_MPI_SEQ_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(2D_MPI_OMP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(2D_MPI_CUDA_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(2D_MPI_HIP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(3D_SEQ_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_3d/++g"),$@/$(dir))

$(3D_OMP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_3d/++g"),$@/$(dir))

$(3D_CUDA_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_3d/++g"),$@/$(dir))

$(3D_HIP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_3d/++g"),$@/$(dir))

$(3D_MPI_SEQ_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_3d/++g"),$@/$(dir))

$(3D_MPI_OMP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_3d/++g"),$@/$(dir))

$(3D_MPI_CUDA_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_3d/++g"),$@/$(dir))

$(3D_MPI_HIP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_3d/++g"),$@/$(dir))

# 2D Common object files
COMMON_2D_OBJ := ins2d_op.o \
	solvers_op/2d/advection_solver_over_int_op.o \
  solvers_op/2d/advection_solver_op.o \
  solvers_op/2d/ls_solver_op.o \
  solvers/2d/ls_utils/ls_reinit_poly.o \
  solvers_op/2d/mp_ins_solver_over_int_op.o \
  solvers_op/2d/ins_solver_base_op.o \
  solvers_op/2d/ins_solver_op.o \
  solvers_op/2d/mp_ins_solver_op.o \
  solvers_op/2d/ins_solver_over_int_op.o \
  solvers_op/2d/ce_solver_over_int_op.o

# 2D CPU only object files
CPU_2D_OBJ := solvers/2d/ls_utils/ls_reinit.o

# TODO 2D CUDA only object files

# TODO 2D HIP only object files

# 2D Non-MPI only object files
SN_2D_OBJ := solvers/2d/ls_utils/kd_tree.o

# 2D MPI only object files
MPI_2D_OBJ := solvers/2d/ls_utils/kd_tree_mpi.o

# 3D Common object files
COMMON_3D_OBJ := ins3d_op.o \
	solvers_op/3d/advection_solver_op.o \
  solvers_op/3d/ins_solver_base_op.o \
  solvers_op/3d/ins_solver_op.o \
  solvers_op/3d/ls_solver_op.o \
  solvers/3d/ls_utils/poly_approx.o \
  solvers/3d/ls_utils/kd_tree.o \
  solvers_op/3d/mp_ins_solver_op.o

# 3D CPU only object files
CPU_3D_OBJ := solvers/3d/ls_solver_cpu.o

# TODO 3D CUDA only object files

# TODO 3D HIP only object files

# 3D MPI only object files
MPI_3D_OBJ := solvers/3d/ls_utils/kd_tree_mpi.o

# Generic rules 2D OpenMP
$(2D_OMP_OBJ_DIR)/%.o: gen_2d/%.cpp | $(2D_OMP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMP_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(INS_INC) -c $< -o $@

# Generic rules 2D MPI + OpenMP
$(2D_MPI_SEQ_OBJ_DIR)/%.o: gen_2d/%.cpp | $(2D_MPI_SEQ_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(SEQ_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 2D MPI + OpenMP
$(2D_MPI_OMP_OBJ_DIR)/%.o: gen_2d/%.cpp | $(2D_MPI_OMP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMP_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 3D OpenMP
$(3D_OMP_OBJ_DIR)/%.o: gen_3d/%.cpp | $(3D_OMP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMP_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(INS_INC) -c $< -o $@

# Generic rules 3D MPI + OpenMP
$(3D_MPI_SEQ_OBJ_DIR)/%.o: gen_3d/%.cpp | $(3D_MPI_SEQ_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(SEQ_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 3D MPI + OpenMP
$(3D_MPI_OMP_OBJ_DIR)/%.o: gen_3d/%.cpp | $(3D_MPI_OMP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMP_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# 2D OpenMP Binary
2D_OMP_OBJ := $(addprefix $(2D_OMP_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(CPU_2D_OBJ) \
	$(SN_2D_OBJ) \
	openmp/ins2d_kernels.o)
$(BIN)/ins2d_openmp: $(2D_OMP_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_openmp -ldgtoolkit $(OP2_OPENMP_LIBS) $(OPENMP_FLAG) $(COMMON_LIBS) -o $@

# 2D MPI Binary
2D_MPI_SEQ_OBJ := $(addprefix $(2D_MPI_SEQ_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(CPU_2D_OBJ) \
	$(MPI_2D_OBJ) \
	seq/ins2d_seqkernels.o)
$(BIN)/ins2d_mpi: $(2D_MPI_SEQ_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_mpi -ldgtoolkit $(OP2_MPI_LIBS) $(PARTITION_LIB) $(COMMON_LIBS) -o $@

# 2D MPI + OpenMP Binary
2D_MPI_OMP_OBJ := $(addprefix $(2D_MPI_OMP_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(CPU_2D_OBJ) \
	$(MPI_2D_OBJ) \
	openmp/ins2d_kernels.o)
$(BIN)/ins2d_mpi_openmp: $(2D_MPI_OMP_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_mpi_openmp -ldgtoolkit $(OP2_MPI_LIBS) $(PARTITION_LIB) $(OPENMP_FLAG) $(COMMON_LIBS) -o $@

# 3D OpenMP Binary
3D_OMP_OBJ := $(addprefix $(3D_OMP_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(CPU_3D_OBJ) \
	openmp/ins3d_kernels.o)
$(BIN)/ins3d_openmp: $(3D_OMP_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_openmp -ldgtoolkit $(OP2_OPENMP_LIBS) $(OPENMP_FLAG) $(COMMON_LIBS) -o $@

# 3D MPI Binary
3D_MPI_SEQ_OBJ := $(addprefix $(3D_MPI_SEQ_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(CPU_3D_OBJ) \
	$(MPI_3D_OBJ) \
	seq/ins3d_seqkernels.o)
$(BIN)/ins3d_mpi: $(3D_MPI_SEQ_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_mpi -ldgtoolkit $(OP2_MPI_LIBS) $(PARTITION_LIB) $(COMMON_LIBS) -o $@

# 3D MPI + OpenMP Binary
3D_MPI_OMP_OBJ := $(addprefix $(3D_MPI_OMP_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(CPU_3D_OBJ) \
	$(MPI_3D_OBJ) \
	openmp/ins3d_kernels.o)
$(BIN)/ins3d_mpi_openmp: $(3D_MPI_OMP_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_mpi_openmp -ldgtoolkit $(OP2_MPI_LIBS) $(PARTITION_LIB) $(OPENMP_FLAG) $(COMMON_LIBS) -o $@