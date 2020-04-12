
# This makefile builds the QuEST library and links QuESTlink.
# While attempting to accomodate as many platforms and compilers,
# unforeseen problems are inevitable: please email tyson.jones@materials.ox.ac.uk 
# about any errors or complications, or raise an issue on Github.
# This makefile is a small change to that created by Tyson Jones for QuEST, 
# which in turn is based off the makefile by Ania Brown for QuEST.

#======================================================================#
#                                                                      #
#      User settings                                                   #
#                                                                      #
#======================================================================#

# operating system, one of {MACOS, LINUX, WINDOWS}
OS = LINUX

# compiler to use, which should support both C and C++, to be wrapped by GPU/MPI compilers
# this is likely to be one of {g++, clang, ic, cl}
COMPILER = g++

# type of above compiler, one of {GNU, INTEL, CLANG, MSVC}, used for setting compiler flags
COMPILER_TYPE = GNU

# only for WINDOWS: whether OS is 32-bit (x86) or 64-bit (x64). Choose {32, 64}
WINDOWS_ARCH = 64

# hardwares to target: 1 means use, 0 means don't use
MULTITHREADED = 0
GPUACCELERATED = 0

# GPU hardware dependent, lookup at https://developer.nvidia.com/cuda-gpus, write without fullstop
GPU_COMPUTE_CAPABILITY = 61

# whether to suppress the below warnings about compiler compatibility
SUPPRESS_WARNING = 0


#======================================================================#
#                                                                      #
#      Constants                                                       #
#                                                                      #
#======================================================================#

# name of the executable to create
EXE = quest_link

# space-separated names (no file type) of all user source files (.c or .cpp) in the root directory
SOURCES = quest_link quest_templates.tm

# path to QuEST library from root directory
QUEST_DIR = QuEST

# path to WSTP libs from root directory 
WSTP_DIR = WSTP

# path to QuESTlink code from root directory
LINK_DIR = Link

# whether to use single, double or quad floating point precision in the state-vector {1,2,4}
PRECISION = 2


#======================================================================#
#                                                                      #
#      Checking user settings                                          #
#                                                                      #
#======================================================================#

# suppresses all non-gcc output, useful for calling scripts
SILENT = 0

# always allow cleaning without errors or warnings
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), veryclean)
ifneq ($(SILENT), 1)

    # check $OS is correct
    ifneq ($(OS), LINUX)
    ifneq ($(OS), MACOS)
    ifneq ($(OS), WINDOWS)
        $(error OS must be LINUX, MACOS or LINUX)
    endif
    endif
    endif

    # check $COMPILER_TYPE is correct
    ifneq ($(COMPILER_TYPE), CLANG)
    ifneq ($(COMPILER_TYPE), GNU)
    ifneq ($(COMPILER_TYPE), INTEL)
    ifneq ($(COMPILER_TYPE), MSVC)
        $(error COMPILER_TYPE must be one of CLANG, GNU, INTEL or MSVC)
    endif
    endif
    endif
	endif

    # distributed GPU not supported
    ifeq ($(DISTRIBUTED), 1)
    ifeq ($(GPUACCELERATED), 1)
        $(error Distributed GPU acceleration not supported)
    endif
    endif

    # GPU doesn't use threading
    ifeq ($(MULTITHREADED), 1)
    ifeq ($(GPUACCELERATED), 1)
        $(warning GPU acceleration makes no use of multithreading. Disabling the latter...)
        override MULTITHREADED = 0
    endif
    endif
	
    # CLANG compilers don't support threading
    ifeq ($(MULTITHREADED), 1)
    ifeq ($(COMPILER_TYPE), CLANG)
        $(warning Clang does not support multithreading. Disabling...)
        override MULTITHREADED = 0
    endif
    endif
	
	# check PRECISION is valid
    ifneq ($(PRECISION), 1)
    ifneq ($(PRECISION), 2)
    ifneq ($(PRECISION), 4)
        $(error PRECISION must be set to 1, 2 or 4)
    endif
    endif
    endif
	
	# GPU does not support quad precision
    ifeq ($(PRECISION), 4)
    ifeq ($(GPUACCELERATED), 1)
        $(warning GPUs do not support quad precision. Setting PRECISION=2...)
        override PRECISION = 2	
    endif
    endif
	
    # NVCC doesn't support new CLANG compilers
    ifeq ($(GPUACCELERATED), 1)
    ifeq ($(COMPILER_TYPE), CLANG)
    ifeq ($(SUPPRESS_WARNING), 0)
        $(info Some versions of Clang are not NVIDIA-GPU compatible. If compilation fails, try Clang 3.7)
    endif
    endif
    endif

    # NVCC doesn't support GNU compilers on OSX
    ifeq ($(GPUACCELERATED), 1)
    ifeq ($(COMPILER_TYPE), GNU)
    ifeq ($(SUPPRESS_WARNING), 0)
        $(info On some platforms (e.g. OSX), NVIDIA-GPUs are not compatible with GNU compilers. If compilation fails, try an alternative compiler, like Clang 3.7)
    endif
    endif
    endif
		
		# Windows users must set WINDOWS_ARCH as {32, 64}
    ifeq ($(OS), WINDOWS)
    ifneq ($(WINDOWS_ARCH), 32)
    ifneq ($(WINDOWS_ARCH), 64)
        $(error When compiling on WINDOWS, WINDOWS_ARCH must be 32 or 64)
    endif
    endif
    endif

# end of allowed cleaning
endif
endif
endif


#======================================================================#
#                                                                      #
#     Compilation                                                      #
#                                                                      #
#======================================================================#


# 
# --- WSTP path
# 

ifeq ($(OS), WINDOWS)
    WSTP_SRC_DIR = $(WSTP_DIR)/Windows
else ifeq ($(OS), MACOS)
    WSTP_SRC_DIR = $(WSTP_DIR)/MacOS
else ifeq ($(OS), LINUX)
    WSTP_SRC_DIR = $(WSTP_DIR)/Linux
endif

#
# --- libraries
#

ifeq ($(OS), MACOS)
    LIBS = -lm -lc++ $(WSTP_SRC_DIR)/MACOS_libWSTPi4.36.a -framework Foundation
else ifeq ($(OS), WINDOWS)
    LIBS = kernel32.lib user32.lib gdi32.lib $(WSTP_SRC_DIR)/windows_wstp$(WINDOWS_ARCH)i4.lib $(WSTP_SRC_DIR)/windows_wstp$(WINDOWS_ARCH)i4m.lib $(WSTP_SRC_DIR)/windows_wstp$(WINDOWS_ARCH)i4s.lib
else ifeq ($(OS), LINUX)
    LIBS = -lm -ldl -lutil -lpthread -luuid -lrt -lstdc++ $(WSTP_SRC_DIR)/linux_libWSTP64i4.a
    ifeq ($(GPUACCELERATED), 0)
        LIBS := -Wl,--no-as-needed $(LIBS)
    endif
endif


#
# --- QuEST source and include paths
#

QUEST_INCLUDE_DIR = ${QUEST_DIR}/include
QUEST_SRC_DIR = ${QUEST_DIR}/src

QUEST_COMMON_DIR = $(QUEST_SRC_DIR)
ifeq ($(GPUACCELERATED), 1)
    QUEST_INNER_DIR = $(QUEST_SRC_DIR)/GPU
else
    QUEST_INNER_DIR = $(QUEST_SRC_DIR)/CPU
endif
QUESTLINK_INCLUDE = -I${QUEST_INCLUDE_DIR} -I$(QUEST_INNER_DIR) -I$(QUEST_COMMON_DIR) -I$(WSTP_SRC_DIR) -I$(LINK_DIR)


#
# --- wrapper compilers
#

CUDA_COMPILER = nvcc



#
# --- compiler flags
#

# note:
#	several flag names depend not just on the compiler type, but also compiler version
#	the user should update these below. For example:
#
#	- GNU C++ compilers of version < 4.7 use -std=c++0x instead of -std=c++11
#	- INTEL compilers of version < ? use -openmp instead of -qopenmp
#	- INTEL compilers of version < ? won't recognise -diad-disable and -cpu-dispatch
#	- CLANG compilers don't support openmp (threading) at all


# threading flag
ifeq ($(MULTITHREADED), 1)
    ifeq ($(COMPILER_TYPE), GNU)
        THREAD_FLAGS = -fopenmp
    else ifeq ($(COMPILER_TYPE), INTEL)
        THREAD_FLAGS = -qopenmp
    else ifeq ($(COMPILER_TYPE), MSVC)
        THREAD_FLAGS = -openmp
    endif
else
    THREAD_FLAGS =
endif

# c
C_CLANG_FLAGS = -O2 -std=c99 -mavx -Wall -DQuEST_PREC=$(PRECISION)
C_GNU_FLAGS = -O2 -std=c99 -mavx -Wall -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS)
C_INTEL_FLAGS = -O2 -std=c99 -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS)
C_MSVC_FLAGS = -O2 -EHs -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS) -nologo -DDWIN$(WINDOWS_ARCH) -D_WINDOWS -Fo$@

# c++
CPP_CLANG_FLAGS = -O2 -std=c++11 -mavx -Wall -DQuEST_PREC=$(PRECISION)
CPP_GNU_FLAGS = -O2 -std=c++11 -mavx -Wall -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS)
CPP_INTEL_FLAGS = -O2 -std=c++11 -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS)
CPP_MSVC_FLAGS = -O2 -EHs -std:c++latest -DQuEST_PREC=$(PRECISION) $(THREAD_FLAGS) -nologo -DDWIN$(WINDOWS_ARCH) -D_WINDOWS -Fo$@

# wrappers
CPP_CUDA_FLAGS = -O2 -arch=compute_$(GPU_COMPUTE_CAPABILITY) -code=sm_$(GPU_COMPUTE_CAPABILITY) -DQuEST_PREC=$(PRECISION) -ccbin $(COMPILER)

# choose c/c++ flags based on compiler type
ifeq ($(COMPILER_TYPE), CLANG)
    C_FLAGS = $(C_CLANG_FLAGS)
    CPP_FLAGS = $(CPP_CLANG_FLAGS)
else ifeq ($(COMPILER_TYPE), GNU)
    C_FLAGS = $(C_GNU_FLAGS)
    CPP_FLAGS = $(CPP_GNU_FLAGS)
else ifeq ($(COMPILER_TYPE), INTEL)
    C_FLAGS = $(C_INTEL_FLAGS)
    CPP_FLAGS = $(CPP_INTEL_FLAGS)
else ifeq ($(COMPILER_TYPE), MSVC)
    C_FLAGS = $(C_MSVC_FLAGS)
    CPP_FLAGS = $(CPP_MSVC_FLAGS)
endif



#
# --- compiler mode and linker flags 
#

ifeq ($(COMPILER_TYPE), MSVC)
    C_MODE = 
    LINKER = link
    ifeq ($(WINDOWS_ARCH), 32)
        ARCH_FLAG = X86
		else
        ARCH_FLAG = X64
		endif
    LINK_FLAGS = -out:$(EXE).exe -SUBSYSTEM:WINDOWS -nologo -MACHINE:$(ARCH_FLAG)
else
    C_MODE = -x c
    LINKER = $(COMPILER)
    LINK_FLAGS = -o $(EXE)
endif



#
# --- targets
#

OBJ = QuEST.o QuEST_validation.o QuEST_common.o QuEST_qasm.o mt19937ar.o
ifeq ($(GPUACCELERATED), 1)
    OBJ += QuEST_gpu.o
else
    OBJ += QuEST_cpu.o QuEST_cpu_local.o
endif
OBJ += $(addsuffix .o, $(SOURCES))


#
# --- rules
#

# notes:
#	- if $SOURCES appear as both c and c++ files, the c files will be compiled
#	- CUDA won't compile .c files ($COMPILER will), only .cpp and .cu

# GPU
ifeq ($(GPUACCELERATED), 1)

  %.o: %.c
	$(COMPILER) $(C_MODE) $(C_FLAGS) $(QUESTLINK_INCLUDE) -c $<
  %.o: $(QUEST_INNER_DIR)/%.c
	$(COMPILER) $(C_MODE) $(C_FLAGS) $(QUESTLINK_INCLUDE) -c $<
  %.o: $(QUEST_COMMON_DIR)/%.c
	$(COMPILER) $(C_MODE) $(C_FLAGS) $(QUESTLINK_INCLUDE) -c $<
  %.o: $(LINK_DIR)/%.c
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $(QUESTLINK_INCLUDE) $<

  %.o: %.cu
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $(QUESTLINK_INCLUDE) $<
  %.o: $(QUEST_INNER_DIR)/%.cu
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $(QUESTLINK_INCLUDE) $<
	
  %.o: %.cpp quest_templates.tm.cpp
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $(QUESTLINK_INCLUDE) $<
  %.o: $(QUEST_INNER_DIR)/%.cpp
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $(QUESTLINK_INCLUDE) $<
  %.o: $(LINK_DIR)/%.cpp
	$(CUDA_COMPILER) -dc $(CPP_CUDA_FLAGS) $(QUESTLINK_INCLUDE) $<

# CPU
else

  %.o: %.c
	$(COMPILER) $(C_MODE) $(C_FLAGS) $(QUESTLINK_INCLUDE) -c $<
  %.o: $(QUEST_INNER_DIR)/%.c
	$(COMPILER) $(C_MODE) $(C_FLAGS) $(QUESTLINK_INCLUDE) -c $<
  %.o: $(QUEST_COMMON_DIR)/%.c
	$(COMPILER) $(C_MODE) $(C_FLAGS) $(QUESTLINK_INCLUDE) -c $<
	
  %.o: %.cpp quest_templates.tm.cpp
	$(COMPILER) $(CPP_FLAGS) $(QUESTLINK_INCLUDE) -c $<
  %.o: $(QUEST_INNER_DIR)/%.cpp
	$(COMPILER) $(CPP_FLAGS) -c $<
  %.o: $(LINK_DIR)/%.cpp
	$(COMPILER) $(CPP_FLAGS) $(QUESTLINK_INCLUDE) -c $<

endif

#
# --- build
#

# CUDA
ifeq ($(GPUACCELERATED), 1)

  all:	$(OBJ)
		$(CUDA_COMPILER) $(CPP_CUDA_FLAGS) $(QUESTLINK_INCLUDE) -o $(EXE) $(OBJ) $(LIBS)

# C
else

  default:	$(EXE)
  $(EXE):	$(OBJ)
			$(LINKER) $(OBJ) $(LIBS) $(LINK_FLAGS)

endif



#
# --- generate C code from MMA templates 
#

ifeq ($(OS), MACOS)
    PREP = MACOS_wsprep
else ifeq ($(OS), LINUX)
    PREP = linux_wsprep
else ifeq ($(OS), WINDOWS)
    PREP = windows_wsprep.exe
endif

quest_templates.tm.cpp:
	$(WSTP_SRC_DIR)/$(PREP) $(LINK_DIR)/quest_templates.tm -o quest_templates.tm.cpp



#
# --- clean
#

# resolve os remove command
ifeq ($(OS), MACOS)
    REM = /bin/rm -f
    EXE_FN = $(EXE)
else ifeq ($(OS), LINUX)
    REM = /bin/rm -f
    EXE_FN = $(EXE)
else ifeq ($(OS), WINDOWS)
    REM = del
    EXE_FN = $(EXE).exe
endif


# define tidy cmds
.PHONY:		tidy clean veryclean
tidy:
			$(REM) *.o
			$(REM) quest_templates.tm.cpp
clean:
			$(REM) *.o $(EXE_FN)
			$(REM) quest_templates.tm.cpp
veryclean:	clean
			$(REM) *.h~ *.c~ makefile~



#
# --- debug
#
	
print-%:
	@echo $*=$($*)

getvalue-%:
	@echo $($*)



