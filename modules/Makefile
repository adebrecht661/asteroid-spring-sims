export OPENGL=1

include Makefile.defs

LIB+=-L/usr/X11/lib -L../src/
LINKER+=-lrebound
OPT+=-Wl,-rpath,./
INCLUDES+= -I../src/ -I../src_spring/

MYSOURCES=../src_spring/matrix_math.cpp \
		  ../src_spring/springs.cpp \
		  ../src_spring/kepcart.cpp \
		  ../src_spring/heat.cpp \
		  ../src_spring/orb.cpp \
		  ../src_spring/input_spring.cpp \
		  ../src_spring/physics.cpp \
		  ../src_spring/shapes.cpp \
		  ../src_spring/stress.cpp \
		  ../src_spring/output_spring.cpp
MYOBJECTS=$(MYSOURCES:.cpp=.o)
MYHEADERS=$(MYSOURCES:.cpp=.h)

MYPROBLEM=$(PROBLEM)/problem.cpp
MYPROBLEMOBJ=$(MYPROBLEM:.cpp=.o)

EXECUTABLE = $(PROBLEM)/rebound_spring

default: $(EXECUTABLE)

$(EXECUTABLE) : $(MYOBJECTS) $(MYPROBLEMOBJ)
	@echo "Compiling $@ ..."
	$(CPP) $(OPT) $(PREDEF) $(INCLUDES) $(LIB) $(MYOBJECTS) $(MYPROBLEMOBJ) -o $(EXECUTABLE) $(LINKER)
	@echo "REBOUND-SPRING compiled successfully."

%.o : %.cpp
	@echo "Compiling source file $< ..."
	$(CPP) -c $(OPT) $(PREDEF) $(INCLUDES) $(LIB) -o $@ $<

clean:
	@echo "Cleaning up shared library librebound.so ..."
	rm -f librebound.so
	$(MAKE) -C ../src/ clean
	@echo "Cleaning up local directory ..."
	rm -rf ../src_spring/*.o
	rm -vf rebound_spring
