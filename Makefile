SHELL := /bin/bash

.PHONY: default
default: all

# Find out the base directory
CURDIR = $(realpath $(PWD) )

OBJDIR = build
$(shell mkdir -p $(OBJDIR))

NAME = AttractServer_mb
BINARY = $(NAME)

SOURCE_DIR = $(CURDIR)/src
SOURCE_DIR_TEST = $(CURDIR)/test

# choose target
# Show list of all available targets for the current config with cfg_help
# Each config should at least provide the two targets RELEASE and DEBUG.
TARGET ?= RELEASE
TEST ?= OFF
CUDA ?= ON

SOURCE_FOLDERS = $(shell find $(SOURCE_DIR) -maxdepth 2 -type d)
SOURCE_TEST_FOLDERS = $(shell find $(SOURCE_DIR_TEST) -maxdepth 2 -type d)

INCLUDES = $(foreach d, $(SOURCE_FOLDERS), -I$d) 

ifeq ($(TEST), OFF)
#	SOURCES_CPP = $(shell find $(SOURCE_DIR) -name "*.cpp" -and -not -path "*emATTRACT*")
	SOURCES_CPP = $(shell find $(SOURCE_DIR) -name "*.cpp")
#	VPATH = $(SOURCE_DIR):$(SOURCE_DIR)/fileIO:$(SOURCE_DIR)/allocator:$(SOURCE_DIR)/cli:$(SOURCE_DIR)/fileIO:$(SOURCE_DIR)/model:$(SOURCE_DIR)/service:$(SOURCE_DIR)/apps:$(SOURCE_DIR)/commons:$(SOURCE_DIR)/server
	VPATH = $(foreach d, $(SOURCE_FOLDERS), $d:)
else ifeq ($(TEST), ON)
	BINARY = "$(NAME)_TEST"
	EXCLUDE = main.cpp
	SOURCES_CPP = $(shell find $(SOURCE_DIR) -name "*.cpp"  ! -name $(EXCLUDE)) \
		$(shell find $(SOURCE_DIR_TEST) -name "*.cpp")	
	LIBS_TEST = -lgtest -lgmock
	VPATH = $(foreach d, $(SOURCE_TEST_FOLDERS), $d:):$(foreach d, $(SOURCE_FOLDERS), $d:)
	INCLUDES += $(foreach d, $(SOURCE_TEST_FOLDERS), -I$d) 
endif

OBJECTS_CPP = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CPP:.cpp=.o)))



# AttractServer 
CXX = g++
ifeq ($(TARGET), RELEASE)
	OFLAGS = -O3 -DNDEBUG
	FXX_OFLAGS = -O2 
else ifeq ($(TARGET), DEBUG) 
	OFLAGS = -O0 -g -Wall -Wextra
	FXX_OFLAGS = -O0 -g -Wall -Wextra
else
	OFLAGS =
	FXX_OFLAGS = 
endif

CXXFLAGS =  $(OFLAGS) -std=c++11 -fmessage-length=0
#INCLUDES = -I$(CURDIR)/src -I$(CURDIR)/src/fileIO 
LDFLAGS =  #-L...
LIBS =  -lpthread -lrt $(LIBS_TEST) -lboost_program_options -lgfortran -lboost_coroutine -lboost_context -lboost_system

# for testing with different boost lib versions. should be disabled by default.
ifeq (1, 0)
#	INCLUDES += -I/home/uwe/dev/boost_1_55_0
#	LDFLAGS += -L/home/uwe/dev/boost_1_55_0/stage/lib
	INCLUDES += -I/home/uwe/dev/boost_1_54_0
	LDFLAGS += -L/home/uwe/dev/boost_1_54_0/stage/lib
endif

ifeq ($(CUDA), ON)
	OFLAGS += -DCUDA
	LDFLAGS += -L$(CUDADIR)/lib64 -Wno-deprecated-gpu-targets
	LIBS += -lcudart -lnvToolsExt
	INCLUDES += -I$(CUDADIR)/include 
	LXX = /usr/local/cuda/bin/nvcc
else
	LXX = ${CXX}
endif



# prepare nvcc settings
CUDA_CXX = /usr/local/cuda/bin/nvcc
ARCHFLAGS = -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50
ARCHFLAGS2 = -gencode arch=compute_30,code=compute_30 -gencode arch=compute_35,code=compute_35 -gencode arch=compute_50,code=compute_50
SOURCES_CU = $(shell find $(SOURCE_DIR) -name "*.cu")

OBJECTS_CU = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CU:.cu=_cu.o)))
ifeq ($(TARGET), RELEASE)
	CUDA_OFLAGS = ${OFLAGS} -D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED
else ifeq ($(TARGET), DEBUG) 
	CUDA_OFLAGS = -O0 -g -G -DCUDA
else
	CUDA_OFLAGS =
endif
CUDA_CXXFLAGS = ${CUDA_OFLAGS} -std=c++11


FXX = gfortran
FXXFLAGS = $(FXX_OFLAGS) -fcray-pointer
SOURCES_F = $(shell find $(CURDIR)/src/ -name "*.f")
OBJECTS_F = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_F:.f=.o)))


OBJECTS = ${OBJECTS_CPP} $(OBJECTS_F)
ifeq ($(CUDA), ON)
	OBJECTS += ${OBJECTS_CU}
endif

# search directories for source files
#VPATH = $(shell find  src/ -type d)
#VPATH += $(shell find test/ -type d)
#VPATH = $(SOURCE_DIR_TEST):$(SOURCE_DIR):$(SOURCE_DIR)/fileIO



$(BINARY): $(OBJECTS)	
	@echo 'Building target: $@'
	@echo 'Invoking: GCC/NVCC C++ Linker'
	$(LXX) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(OBJECTS_CPP): $(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d
	@echo 'Building file: "$<"'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '

$(OBJDIR)/%.d:	;

$(OBJECTS_F): $(OBJDIR)/%.o: %.f
	@echo 'Building file: "$<"'
	@echo 'Invoking: GCC Fortran Compiler'
	$(FXX) $(FXXFLAGS) -c -o "$@" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '	
	

-include $(OBJECTS_CPP:.o=.d)

$(OBJECTS_CU): $(OBJDIR)/%_cu.o: %.cu 
	@echo 'Building file: $<"'
	@echo 'Invoking: NVCC Compiler'
	$(CUDA_CXX) $(CUDA_CXXFLAGS) $(ARCHFLAGS) $(ARCHFLAGS2) $(INCLUDES) --compile --relocatable-device-code=true -x cu -o "$@" "$<"
	@echo 'Finished building: $<"'
	@echo ' '	

.PHONY: all clean cleanall clean_deps

all: $(BINARY)

clean:
	rm -r $(OBJDIR)
	
cleanall: clean
	find . -type f -name "$(BINARY)*" -delete

clean_deps:
	find $(OBJDIR) -type f -name "*.d" -delete

