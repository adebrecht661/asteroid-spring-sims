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

default: librebound $(PROBLEM)

all: librebound $(PROBLEMS)

$(PROBLEMS): librebound 
	cd modules && make PROBLEM=$@ INC=../$(INC)

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C src/

clean:
	cd modules; make clean
	for PROBLEM in $(PROBLEM); do  \
	  rm -f modules/$$PROBLEM/problem.o; \
	  rm -f modules/$$PROBLEM/rebound_spring; \
	done

allclean: 
	cd src; make clean
	for PROBLEM in $(PROBLEMS); do  \
	  rm -f modules/$$PROBLEM/problem.o; \
	  rm -f modules/$$PROBLEM/rebound_spring; \
	done