include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

TopOpt: TopOpt.o
	-${CLINKER} -o TopOpt TopOpt.o ${PETSC_LIB}
	${RM} TopOpt.o
TopOpt2: TopOpt2.o
	-${CLINKER} -o TopOpt2 TopOpt2.o ${PETSC_LIB}
	${RM} TopOpt2.o
