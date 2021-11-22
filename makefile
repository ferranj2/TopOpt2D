include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test
CFLAGS += -pedantic -std=c99

TopOpt: TopOpt.o
	-${CLINKER} -o TopOpt TopOpt.o ${PETSC_LIB}
	${RM} TopOpt.o
TopOpt2: TopOpt2.o
	-${CLINKER} -o TopOpt2 TopOpt2.o ${PETSC_SYS_LIB}
	${RM} TopOpt2.o
TopOpt2filter: TopOpt2filter.o
	-${CLINKER} -o TopOpt2filter TopOpt2filter.o ${PETSC_LIB}
	${RM} TopOpt2filter.o
TestMatGenAndPlot: TestMatGenAndPlot.o  chkopts
        -${CLINKER} -o TestMatGenAndPlot TestMatGenAndPlot.o  ${PETSC_SYS_LIB}
        ${RM} TestMatGenAndPlot.o
CleanImage:
	rm -f tmpDat/*.dat rm -f tmpAni/*.jpg

