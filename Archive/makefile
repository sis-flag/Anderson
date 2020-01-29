SLEPC_DIR=/Users/lihuiyuan/GNU/slepc-3.7.3
PETSC_DIR=/Users/lihuiyuan/GNU/petsc
PETSC_ARCH=arch-darwin-c-debug
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common


test1d: schrand.o
	-${CLINKER} -o test1d schrand.o  ${SLEPC_EPS_LIB}
	${RM} schrand.o
	${MPIEXEC} -n 2 ./test1d -m 2 -n 100
#	${RM} test1d
