MAKEFILES_DIR != dirname $(realpath \
  $(word $(words $(MAKEFILE_LIST)), $(MAKEFILE_LIST)))
ROOT_DIR != realpath $(MAKEFILES_DIR)

INC := -I$(ROOT_DIR)/include
BIN := $(ROOT_DIR)/bin
OBJ := $(ROOT_DIR)/obj
CODE_GEN_DIR := $(ROOT_DIR)/code_gen
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

COMMON_LIBS = $(OPENBLAS_LIB) $(PETSC_LIB) $(ARMA_LIB) $(HDF5_LIB)

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
HIP_COMPILE_DEFS := -DOP2_DG_HIP -DDG_OP2_SOA -DINS_HIP
CUDA_COMPILE_DEFS := -DOP2_DG_CUDA -DDG_OP2_SOA -DINS_CUDA
MPI_COMPILER_DEFS := -DDG_MPI -DINS_MPI
INS_INC := $(ARMA_INC) $(INC) $(OP2_INC) $(OPENBLAS_INC) $(PETSC_INC) $(OP2_DG_TOOLKIT_INC) $(INIPP_INC) $(CDT_INC)

ifeq ($(BUILD_WITH_HYPRE),1)
	INS_INC := $(INS_INC) $(HYPRE_INC)
	COMMON_COMPILE_DEFS_2D := $(COMMON_COMPILE_DEFS_2D) -DINS_BUILD_WITH_HYPRE
	COMMON_COMPILE_DEFS_3D := $(COMMON_COMPILE_DEFS_3D) -DINS_BUILD_WITH_HYPRE
	COMMON_LIBS := $(COMMON_LIBS) $(HYPRE_LIB)
endif

all: cpu_all cuda_all

cpu_all: cpu_2d_all cpu_3d_all
cpu_2d_all: cpu_omp_2d cpu_mpi_seq_2d cpu_mpi_omp_2d
cpu_3d_all: cpu_omp_3d cpu_mpi_seq_3d cpu_mpi_omp_3d

cpu_omp_2d: $(BIN)/ins2d_openmp
cpu_mpi_seq_2d: $(BIN)/ins2d_mpi
cpu_mpi_omp_2d: $(BIN)/ins2d_mpi_openmp

cpu_omp_3d: $(BIN)/ins3d_openmp
cpu_mpi_seq_3d: $(BIN)/ins3d_mpi
cpu_mpi_omp_3d: $(BIN)/ins3d_mpi_openmp

cuda_all: cuda_2d_all cuda_3d_all
cuda_2d_all: cuda_2d mpi_cuda_2d
cuda_3d_all: cuda_3d mpi_cuda_3d

cuda_2d: $(BIN)/ins2d_cuda
mpi_cuda_2d: $(BIN)/ins2d_mpi_cuda

cuda_3d: $(BIN)/ins3d_cuda
mpi_cuda_3d: $(BIN)/ins3d_mpi_cuda

hip_all: hip_2d_all hip_3d_all
hip_2d_all: hip_2d mpi_hip_2d
hip_3d_all: hip_3d mpi_hip_3d

hip_2d: $(BIN)/ins2d_hip
mpi_hip_2d: $(BIN)/ins2d_mpi_hip

hip_3d: $(BIN)/ins3d_hip
mpi_hip_3d: $(BIN)/ins3d_mpi_hip

codegen: $(CODE_GEN_DIR)

clean:
	-rm -rf $(OBJ)
	-rm -rf $(BIN)
	-rm -rf $(CODE_GEN_DIR)

$(CODE_GEN_DIR):
	OP2_TRANSLATOR=$(OP2_TRANSLATOR) ORDER=$(ORDER) SOA=$(SOA) ./code_gen.sh

$(BIN):
	@mkdir -p $@

$(OBJ):
	@mkdir -p $@

$(2D_SEQ_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+gen_2d/++g"),$@/$(dir))

$(2D_OMP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_2d/++g"),$@/$(dir))

$(2D_CUDA_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_2d/++g"),$@/$(dir))

$(2D_HIP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_2d/++g"),$@/$(dir))

$(2D_MPI_SEQ_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_2d/++g"),$@/$(dir))

$(2D_MPI_OMP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_2d/++g"),$@/$(dir))

$(2D_MPI_CUDA_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_2d/++g"),$@/$(dir))

$(2D_MPI_HIP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_2d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_2d/++g"),$@/$(dir))

$(3D_SEQ_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_3d/++g"),$@/$(dir))

$(3D_OMP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_3d/++g"),$@/$(dir))

$(3D_CUDA_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_3d/++g"),$@/$(dir))

$(3D_HIP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_3d/++g"),$@/$(dir))

$(3D_MPI_SEQ_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_3d/++g"),$@/$(dir))

$(3D_MPI_OMP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_3d/++g"),$@/$(dir))

$(3D_MPI_CUDA_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_3d/++g"),$@/$(dir))

$(3D_MPI_HIP_OBJ_DIR): $(OBJ)
	@mkdir -p $@
	@mkdir -p $(foreach dir,$(shell find $(CODE_GEN_DIR)/gen_3d -maxdepth 3 -mindepth 1 | grep -v "\." | sed "s+$(CODE_GEN_DIR)/gen_3d/++g"),$@/$(dir))

# 2D Common object files
COMMON_2D_OBJ := ins2d_op.o \
  solvers_op/2d/advection_solver_op.o \
	solvers_op/2d/diffusion_solver_op.o \
	solvers_op/2d/advec_diff_solver_op.o \
  solvers_op/2d/ls_solver_op.o \
  solvers/2d/ls_utils/ls_reinit_poly.o \
	solvers/2d/ls_utils/kd_tree.o \
  solvers_op/2d/ins_solver_base_op.o \
  solvers_op/2d/ins_solver_op.o \
	solvers_op/2d/ins_temperature_solver_op.o \
  solvers_op/2d/mp_ins_solver_op.o \
	solvers_op/2d/ls_driver_op.o \
	measurements/2d/measurement.o \
	measurements_op/2d/lift_drag_cylinder_op.o \
	measurements_op/2d/l2_error_vortex_op.o \
	measurements_op/2d/min_max_interface_op.o \
	measurements_op/2d/min_max_pressure_op.o \
	measurements_op/2d/ls_advec_error_op.o

# 2D CPU only object files
CPU_2D_OBJ := solvers/2d/ls_utils/ls_reinit.o

# 2D CUDA only object files
CUDA_2D_OBJ := solvers/2d/ls_utils/ls_reinit_gpu.o

# 2D HIP only object files
HIP_2D_OBJ := solvers/2d/ls_utils/ls_reinit_hip.o

# 2D MPI only object files
MPI_2D_OBJ := solvers/2d/ls_utils/kd_tree_mpi.o

# 3D Common object files
COMMON_3D_OBJ := ins3d_op.o \
	solvers_op/3d/advection_solver_op.o \
	solvers_op/3d/diffusion_solver_op.o \
  solvers_op/3d/ins_solver_base_op.o \
  solvers_op/3d/ins_solver_op.o \
  solvers_op/3d/ls_solver_op.o \
  solvers/3d/ls_utils/poly_approx.o \
  solvers/3d/ls_utils/kd_tree.o \
  solvers_op/3d/mp_ins_solver_op.o \
	solvers_op/3d/ls_driver_op.o \
	measurements/3d/measurement.o \
	measurements_op/3d/enstrophy_op.o

# 3D CPU only object files
CPU_3D_OBJ := solvers/3d/ls_solver_cpu.o

# 3D CUDA only object files
CUDA_3D_OBJ := solvers/3d/ls_solver_gpu.o

# 3D HIP only object files
HIP_3D_OBJ := solvers/3d/ls_solver_hip.o

# 3D MPI only object files
MPI_3D_OBJ := solvers/3d/ls_utils/kd_tree_mpi.o

# Generic rules 2D OpenMP
$(2D_OMP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cpp | $(2D_OMP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMP_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(INS_INC) -c $< -o $@

# Generic rules 2D CUDA
$(2D_CUDA_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cpp | $(2D_CUDA_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CUDA_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(CUDA_COMPILE_DEFS) $(INS_INC) -c $< -o $@

$(2D_CUDA_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cu | $(2D_CUDA_OBJ_DIR)
	$(NVCC) $(NVCC_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(CUDA_COMPILE_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 2D HIP
$(2D_HIP_OBJ_DIR)/hip/ins2d_kernels.o: $(CODE_GEN_DIR)/gen_2d/hip/ins2d_kernels.cpp | $(2D_HIP_OBJ_DIR)
	$(HIPCC) $(HIP_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(HIP_COMPILE_DEFS) $(INS_INC) -c $< -o $@

$(2D_HIP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cpp | $(2D_HIP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(HIP_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(HIP_COMPILE_DEFS) $(INS_INC) -c $< -o $@

$(2D_HIP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.hip | $(2D_HIP_OBJ_DIR)
	$(HIPCC) $(HIP_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(HIP_COMPILE_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 2D MPI + OpenMP
$(2D_MPI_SEQ_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cpp | $(2D_MPI_SEQ_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(SEQ_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 2D MPI + OpenMP
$(2D_MPI_OMP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cpp | $(2D_MPI_OMP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMP_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 2D MPI + CUDA
$(2D_MPI_CUDA_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cpp | $(2D_MPI_CUDA_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CUDA_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(CUDA_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

$(2D_MPI_CUDA_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cu | $(2D_MPI_CUDA_OBJ_DIR)
	$(NVCC) $(NVCC_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(CUDA_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 2D MPI + HIP
$(2D_MPI_HIP_OBJ_DIR)/hip/ins2d_kernels.o: $(CODE_GEN_DIR)/gen_2d/hip/ins2d_kernels.cpp | $(2D_MPI_HIP_OBJ_DIR)
	$(HIPCC) $(HIP_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(HIP_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

$(2D_MPI_HIP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cpp | $(2D_MPI_HIP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(HIP_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(HIP_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

$(2D_MPI_HIP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_2d/%.cu | $(2D_MPI_HIP_OBJ_DIR)
	$(HIPCC) $(HIP_FLAGS) $(COMMON_COMPILE_DEFS_2D) $(HIP_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 3D OpenMP
$(3D_OMP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cpp | $(3D_OMP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMP_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(INS_INC) -c $< -o $@

# Generic rules 3D CUDA
$(3D_CUDA_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cpp | $(3D_CUDA_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CUDA_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(CUDA_COMPILE_DEFS) $(INS_INC) -c $< -o $@

$(3D_CUDA_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cu | $(3D_CUDA_OBJ_DIR)
	$(NVCC) $(NVCC_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(CUDA_COMPILE_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 3D HIP
$(3D_HIP_OBJ_DIR)/hip/ins3d_kernels.o: $(CODE_GEN_DIR)/gen_3d/hip/ins3d_kernels.cpp | $(3D_HIP_OBJ_DIR)
	$(HIPCC) $(HIP_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(HIP_COMPILE_DEFS) $(INS_INC) -c $< -o $@

$(3D_HIP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cpp | $(3D_HIP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(HIP_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(HIP_COMPILE_DEFS) $(INS_INC) -c $< -o $@

$(3D_HIP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.hip | $(3D_HIP_OBJ_DIR)
	$(HIPCC) $(HIP_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(HIP_COMPILE_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 3D MPI + OpenMP
$(3D_MPI_SEQ_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cpp | $(3D_MPI_SEQ_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(SEQ_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 3D MPI + OpenMP
$(3D_MPI_OMP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cpp | $(3D_MPI_OMP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMP_CPU_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 3D MPI + CUDA
$(3D_MPI_CUDA_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cpp | $(3D_MPI_CUDA_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CUDA_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(CUDA_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

$(3D_MPI_CUDA_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cu | $(3D_MPI_CUDA_OBJ_DIR)
	$(NVCC) $(NVCC_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(CUDA_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# Generic rules 3D MPI + HIP
$(3D_MPI_HIP_OBJ_DIR)/hip/ins3d_kernels.o: $(CODE_GEN_DIR)/gen_3d/hip/ins3d_kernels.cpp | $(3D_MPI_HIP_OBJ_DIR)
	$(HIPCC) $(HIP_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(HIP_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

$(3D_MPI_HIP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cpp | $(3D_MPI_HIP_OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(HIP_COMPILER_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(HIP_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

$(3D_MPI_HIP_OBJ_DIR)/%.o: $(CODE_GEN_DIR)/gen_3d/%.cu | $(3D_MPI_HIP_OBJ_DIR)
	$(HIPCC) $(HIP_FLAGS) $(COMMON_COMPILE_DEFS_3D) $(HIP_COMPILE_DEFS) $(MPI_COMPILER_DEFS) $(INS_INC) -c $< -o $@

# 2D OpenMP Binary
2D_OMP_OBJ := $(addprefix $(2D_OMP_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(CPU_2D_OBJ) \
	openmp/ins2d_kernels.o)
$(BIN)/ins2d_openmp: $(2D_OMP_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_openmp -ldgtoolkit $(OP2_OPENMP_LIBS) $(OPENMP_FLAG) $(COMMON_LIBS) $(EXTRA_LIBS_CPU) -o $@

# 2D CUDA Binary
2D_CUDA_OBJ := $(addprefix $(2D_CUDA_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(CUDA_2D_OBJ) \
	cuda/ins2d_kernels.o)
$(BIN)/ins2d_cuda: $(2D_CUDA_OBJ) | $(BIN)
	$(NVCC) $(NVCC_FLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_cuda -ldgtoolkit $(OP2_CUDA_LIBS) $(COMMON_LIBS) $(MPI_LIB) $(CUBLAS_LIB) $(EXTRA_LIBS_CUDA) -o $@

# 2D HIP Binary
2D_HIP_OBJ := $(addprefix $(2D_HIP_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(HIP_2D_OBJ) \
	hip/ins2d_kernels.o)
$(BIN)/ins2d_hip: $(2D_HIP_OBJ) | $(BIN)
	$(HIPCC) $(HIP_FLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_hip -ldgtoolkit $(OP2_HIP_LIBS) $(COMMON_LIBS) $(MPI_LIB) $(HIPBLAS_LIB) $(EXTRA_LIBS_HIP) -o $@

# 2D MPI Binary
2D_MPI_SEQ_OBJ := $(addprefix $(2D_MPI_SEQ_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(CPU_2D_OBJ) \
	$(MPI_2D_OBJ) \
	seq/ins2d_seqkernels.o)
$(BIN)/ins2d_mpi: $(2D_MPI_SEQ_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_mpi -ldgtoolkit $(OP2_MPI_LIBS) $(PARTITION_LIB) $(COMMON_LIBS) $(EXTRA_LIBS_CPU) -o $@

# 2D MPI + OpenMP Binary
2D_MPI_OMP_OBJ := $(addprefix $(2D_MPI_OMP_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(CPU_2D_OBJ) \
	$(MPI_2D_OBJ) \
	openmp/ins2d_kernels.o)
$(BIN)/ins2d_mpi_openmp: $(2D_MPI_OMP_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_mpi_openmp -ldgtoolkit $(OP2_MPI_LIBS) $(PARTITION_LIB) $(OPENMP_FLAG) $(COMMON_LIBS) $(EXTRA_LIBS_CPU) -o $@

# 2D MPI + CUDA Binary
2D_MPI_CUDA_OBJ := $(addprefix $(2D_MPI_CUDA_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(CUDA_2D_OBJ) \
	$(MPI_2D_OBJ) \
	cuda/ins2d_kernels.o)
$(BIN)/ins2d_mpi_cuda: $(2D_MPI_CUDA_OBJ) | $(BIN)
	$(NVCC) $(NVCC_FLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_mpi_cuda -ldgtoolkit $(OP2_MPI_CUDA_LIBS) $(PARTITION_LIB) $(COMMON_LIBS) $(MPI_LIB) $(CUBLAS_LIB) $(EXTRA_LIBS_CUDA) -o $@

# 2D MPI + HIP Binary
2D_MPI_HIP_OBJ := $(addprefix $(2D_MPI_HIP_OBJ_DIR)/,\
	$(COMMON_2D_OBJ) \
	$(HIP_2D_OBJ) \
	$(MPI_2D_OBJ) \
	hip/ins2d_kernels.o)
$(BIN)/ins2d_mpi_hip: $(2D_MPI_HIP_OBJ) | $(BIN)
	$(HIPCC) $(HIP_FLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_2d_mpi_hip -ldgtoolkit $(OP2_MPI_HIP_LIBS) $(PARTITION_LIB) $(COMMON_LIBS) $(MPI_LIB) $(HIPBLAS_LIB) $(EXTRA_LIBS_HIP) -o $@

# 3D OpenMP Binary
3D_OMP_OBJ := $(addprefix $(3D_OMP_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(CPU_3D_OBJ) \
	openmp/ins3d_kernels.o)
$(BIN)/ins3d_openmp: $(3D_OMP_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_openmp -ldgtoolkit $(OP2_OPENMP_LIBS) $(OPENMP_FLAG) $(COMMON_LIBS) $(EXTRA_LIBS_CPU) -o $@

# 3D CUDA Binary
3D_CUDA_OBJ := $(addprefix $(3D_CUDA_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(CUDA_3D_OBJ) \
	cuda/ins3d_kernels.o)
$(BIN)/ins3d_cuda: $(3D_CUDA_OBJ) | $(BIN)
	$(NVCC) $(NVCC_FLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_cuda -ldgtoolkit $(OP2_CUDA_LIBS) $(COMMON_LIBS) $(MPI_LIB) $(CUBLAS_LIB) $(EXTRA_LIBS_CUDA) -o $@

# 3D HIP Binary
3D_HIP_OBJ := $(addprefix $(3D_HIP_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(HIP_3D_OBJ) \
	hip/ins3d_kernels.o)
$(BIN)/ins3d_hip: $(3D_HIP_OBJ) | $(BIN)
	$(HIPCC) $(HIP_FLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_hip -ldgtoolkit $(OP2_HIP_LIBS) $(COMMON_LIBS) $(MPI_LIB) $(HIPBLAS_LIB) $(EXTRA_LIBS_HIP) -o $@

# 3D MPI Binary
3D_MPI_SEQ_OBJ := $(addprefix $(3D_MPI_SEQ_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(CPU_3D_OBJ) \
	$(MPI_3D_OBJ) \
	seq/ins3d_seqkernels.o)
$(BIN)/ins3d_mpi: $(3D_MPI_SEQ_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_mpi -ldgtoolkit $(OP2_MPI_LIBS) $(PARTITION_LIB) $(COMMON_LIBS) $(EXTRA_LIBS_CPU) -o $@

# 3D MPI + OpenMP Binary
3D_MPI_OMP_OBJ := $(addprefix $(3D_MPI_OMP_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(CPU_3D_OBJ) \
	$(MPI_3D_OBJ) \
	openmp/ins3d_kernels.o)
$(BIN)/ins3d_mpi_openmp: $(3D_MPI_OMP_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_mpi_openmp -ldgtoolkit $(OP2_MPI_LIBS) $(PARTITION_LIB) $(OPENMP_FLAG) $(COMMON_LIBS) $(EXTRA_LIBS_CPU) -o $@

# 3D MPI + CUDA Binary
3D_MPI_CUDA_OBJ := $(addprefix $(3D_MPI_CUDA_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(CUDA_3D_OBJ) \
	$(MPI_3D_OBJ) \
	cuda/ins3d_kernels.o)
$(BIN)/ins3d_mpi_cuda: $(3D_MPI_CUDA_OBJ) | $(BIN)
	$(NVCC) $(NVCC_FLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_mpi_cuda -ldgtoolkit $(OP2_MPI_CUDA_LIBS) $(PARTITION_LIB) $(COMMON_LIBS) $(MPI_LIB) $(CUBLAS_LIB) $(EXTRA_LIBS_CUDA) -o $@

# 3D MPI + HIP Binary
3D_MPI_HIP_OBJ := $(addprefix $(3D_MPI_HIP_OBJ_DIR)/,\
	$(COMMON_3D_OBJ) \
	$(HIP_3D_OBJ) \
	$(MPI_3D_OBJ) \
	hip/ins3d_kernels.o)
$(BIN)/ins3d_mpi_hip: $(3D_MPI_HIP_OBJ) | $(BIN)
	$(HIPCC) $(HIP_FLAGS) $^ -L$(OP2_DG_TOOLKIT_DIR)/lib -lop2dgtoolkit_3d_mpi_hip -ldgtoolkit $(OP2_MPI_HIP_LIBS) $(PARTITION_LIB) $(COMMON_LIBS) $(MPI_LIB) $(HIPBLAS_LIB) $(EXTRA_LIBS_HIP) -o $@
