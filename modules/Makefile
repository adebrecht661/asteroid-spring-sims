export OPENGL=1

include Makefile.defs

LIB+=-L../src/ -L../libconfig/lib/.libs/ -L../../libconfig/lib/.libs/
LINKER+=-lrebound -lconfig++
INCLUDES+= -I../src/ -I../src_spring/ -I../libconfig/lib/ -I../../libconfig/lib/
LDFLAGS=-Wl,-R -Wl,../../src/ -Wl,-R -Wl,../../libconfig/lib/.libs/

MYSOURCES=../src_spring/matrix_math.cpp \
		  ../src_spring/springs.cpp \
		  ../src_spring/kepcart.cpp \
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

tests: $(MYOBJECTS) ../src_spring/tests.cpp
	$(CPP) -o tests/tests ../src_spring/tests.cpp $(LIB) $(INCLUDES) $(LDFLAGS) -I../googletest/googletest/include/ ../googletest/lib/libgtest.a $(MYOBJECTS) -pthread $(LINKER)

$(EXECUTABLE) : $(MYOBJECTS) $(MYPROBLEMOBJ)
	@echo "Compiling $@ ..."
	$(CPP) $(OPT) $(PREDEF) $(INCLUDES) $(LIB) $(LDFLAGS) $(MYOBJECTS) $(MYPROBLEMOBJ) -o $(EXECUTABLE) $(LINKER)
	@echo "REBOUND-SPRING compiled successfully."

%.o : %.cpp
	@echo "Compiling source file $< ..."
	$(CPP) -c $(OPT) $(PREDEF) $(INCLUDES) $(LIB) $(LDFLAGS) -o $@ $< $(LINKER)
