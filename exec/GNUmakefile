EXEC = test
SRCDIR = ../src
BUILDDIR = ../build

FC = ifort
FCFLAGS = -g -O0 -CB

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(BUILDDIR)/test.o $(BUILDDIR)/MGPoisson2D.o
	$(FC) $(FCFLAGS) $(BUILDDIR)/test.o $(BUILDDIR)/MGPoisson2D.o -o $(EXEC)

$(BUILDDIR)/test.o: test.f90 $(BUILDDIR)/MGPoisson2D.o
	$(FC) $(FCFLAGS) -c test.f90 -o $(BUILDDIR)/test.o -I$(BUILDDIR)

$(BUILDDIR)/MGPoisson2D.o: $(SRCDIR)/MGPoisson2D.f90
	$(FC) $(FCFLAGS) -c $(SRCDIR)/MGPoisson2D.f90 -o $(BUILDDIR)/MGPoisson2D.o -module $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR)/* $(EXEC)
