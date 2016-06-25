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
	VPATH = $(SOURCE_DIR):$(SOURCE_DIR_TEST)
else ifeq ($(TEST), ON)
	BINARY = "$(NAME)_TEST"
	EXCLUDE = main.cpp
	SOURCES_CPP = $(shell find $(SOURCE_DIR) -name "*.cpp" -and -not -name "$(EXCLUDE)") \
		$(shell find $(SOURCE_DIR_TEST) -name "*.cpp")
	LIBS_TEST = -lgtest -lgmock
	VPATH = $(SOURCE_DIR_TEST):$(SOURCE_DIR)
endif

OBJECTS_CPP = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CPP:.cpp=.o)))
OBJECTS = $(OBJECTS_CPP)

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
INCLUDES = -I$(CUDADIR)/include -I$(CURDIR)/src
LDFLAGS = #-L...
LIBS = -lpthread $(LIBS_TEST)

# search directories for source files
#VPATH = $(shell find ../src/ -type d)
#VPATH = $(SOURCE_DIR_TEST):$(SOURCE_DIR)

$(BINARY): $(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(OBJECTS): $(OBJDIR)/%.o: %.cpp $(OBJDIR)/%.d
	@echo 'Building file: "$<"'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '

$(OBJDIR)/%.d: ;	
	
-include $(OBJECTS:.o=.d)
	
.PHONY: all clean cleanall clean_deps

all: $(BINARY)

clean:
	rm -r $(OBJDIR)
	
cleanall: clean
	find . -type f -name "$(BINARY)*" -delete

clean_deps:
	find $(OBJDIR) -type f -name "*.d" -delete

