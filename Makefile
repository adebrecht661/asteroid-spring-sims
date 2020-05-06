# Makefile for spring code
# The behavior of make and make clean is dependent on whether there is a src/modules/Problem directory.
# 	If there is, than it will only build (or clean) that problem module directory
# 	If there is not, than it will build or clean every subdirectory in src/modules
# You can explicitly set the Problem module directory by passing PROBLEM=<dir> to make or make clean
# You can also build a particular problem module by running make <dir>, but to clean it you will need to 

# Get list of problem directories
PROBLEMS := $(patsubst modules/,,$(patsubst modules/%/,%,$(sort $(dir $(wildcard modules/*/)))))

# Check for whether user specified the value of Problem using make PROBLEM=<dir>
ifndef PROBLEM

# Check if there is a /src/modules/Problem directory
  ifeq (Problem,$(filter Problem,$(PROBLEMS)))
    PROBLEM=Problem  #build Problem
  else
    PROBLEM=$(PROBLEMS) #build everything
  endif
endif

NUM_PROC ?= 1
log_file ?= $(shell date -I).log

default: rebound config $(PROBLEM)

all: rebound config $(PROBLEMS)
 
tests: rebound config gtest
	cd modules && make tests

$(PROBLEMS): rebound config
	cd modules && make PROBLEM=$@ INC=../$(INC)

rebound: 
	@echo "Compiling shared library librebound ..."
	$(MAKE) -C src/
	
config:
	@echo "Compiling shared library libconfig ..."
	cd libconfig; if [ ! -f "config.status" ] ; then ./configure; fi; make
	
gtest:
	@echo "Compiling static library gtest ..."
	cd googletest && cmake CMakeLists.txt && make

clean:
	@echo "Removing spring objects ..." 
	cd src_spring; rm -f *.o
	@echo "Removing spring problems ..."
	for PROBLEM in $(PROBLEM); do  \
	  rm -f modules/$$PROBLEM/problem.o; \
	  rm -f modules/$$PROBLEM/rebound_spring; \
	done

allclean:
	cd libconfig; make clean
	cd src; make clean
	@echo "Removing spring objects ..."
	cd src_spring; rm -f *.o
	@echo "Removing spring problems ..."
	for PROBLEM in $(PROBLEM); do  \
	  rm -f modules/$$PROBLEM/problem.o; \
	  rm -f modules/$$PROBLEM/rebound_spring; \
	  rm -f modules/$$PROBLEM/librebound.so; \
	done