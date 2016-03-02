ALL: hello debug vector poisson2d poisson3d
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/variables

hello : hello.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}
	rm -f hello.o

debug : debug.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}
	rm -f debug.o

vector : vector.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}
	rm -f vector.o

poisson2d : poisson2d.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}
	rm -f poisson2d.o

poisson3d : poisson3d.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}
	rm -f poisson3d.o
