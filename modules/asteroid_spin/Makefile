export OPENGL=1

include Makefile.defs

LIB+=-L/usr/X11/lib
OPT+=-I/usr/X11/include

MYSOURCES=../../src_spring/matrix_math.cpp \
		  ../../src_spring/spring.cpp \
		  ../../src_spring/kepcart.cpp \
		  ../../src_spring/heat.cpp \
		  ../../src_spring/orb.cpp \
		  ../../src_spring/m_output.cpp
MYOBJECTS=$(SOURCES:.c=.o)
MYHEADERS=$(SOURCES:.c=.h)

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../../src/ -I./ -Wl,-rpath,./ $(OPT) $(PREDEF) -I../../src_spring/ $(MYSOURCES) problem.cpp -L. -lrebound $(LIB) -o rebound_spring
	@echo ""
	@echo "REBOUND-SPRING compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../../src/ 
	@-rm -f librebound.so
	@ln -s ../../src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -rf ../../src_spring/*.o
	@-rm -vf rebound_spring
