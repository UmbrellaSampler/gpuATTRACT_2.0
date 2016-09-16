SHELL := /bin/bash

.PHONY: default
default: all

# Find out the base directory
CURDIR = $(realpath $(PWD) )

OBJDIR = build
$(shell mkdir -p $(OBJDIR))

NAME = AttractServer
BINARY = $(NAME)

SOURCE_DIR = $(CURDIR)/src
SOURCE_DIR_TEST = $(CURDIR)/test

# choose target
# Show list of all available targets for the current config with cfg_help
# Each config should at least provide the two targets RELEASE and DEBUG.
TARGET ?= RELEASE
TEST ?= OFF

LIBS_TEST = 
ifeq ($(TEST), OFF)
	SOURCES_CPP = $(shell find $(SOURCE_DIR) -name "*.cpp")
	VPATH = $(SOURCE_DIR):$(SOURCE_DIR)/fileIO
else ifeq ($(TEST), ON)
	BINARY = "$(NAME)_TEST"
	EXCLUDE = main.cpp
	SOURCES_CPP = $(shell find $(SOURCE_DIR) -name "*.cpp" -and -not -name "$(EXCLUDE)") \
		$(shell find $(SOURCE_DIR_TEST) -name "*.cpp")
	LIBS_TEST = -lgtest -lgmock
	VPATH = $(SOURCE_DIR_TEST):$(SOURCE_DIR):$(SOURCE_DIR)/fileIO
endif

OBJECTS_CPP = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CPP:.cpp=.o)))



# AttractServer 
CXX = g++
ifeq ($(TARGET), RELEASE)
	OFLAGS = -O3 -DNDEBUG
else ifeq ($(TARGET), DEBUG) 
	OFLAGS = -O0 -g -Wall -Wextra
else
	OFLAGS =
endif

CXXFLAGS =  $(OFLAGS) -std=c++11 -fmessage-length=0
INCLUDES = -I$(CUDADIR)/include -I$(CURDIR)/src -I$(CURDIR)/src/fileIO 
LDFLAGS = #-L...
LIBS = -lpthread -lrt $(LIBS_TEST) -lboost_program_options


ifeq ($(CUDA), ON)
	OFLAGS += -DCUDA
	LDFLAGS += -L$(CUDADIR)/lib64
	LIBS += -lcudart
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

OBJECTS = ${OBJECTS_CPP}
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

