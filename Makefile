# librairies de SuiteSparse
L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a 
L3 = SuiteSparse/AMD/Lib/libamd.a 
L4 = SuiteSparse/CAMD/Lib/libcamd.a  
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a 
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a 
L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a
LIB = $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack
# M1 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_intelmpi_lp64.a
# M2 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_intelmpi_ilp64.so
# M3 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_openmpi_ilp64.so
# M4 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_openmpi_lp64.so
# M5 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_sgimpt_ilp64.so
# M6 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_sgimpt_lp64.so 
# M7 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blas95_ilp64.a
# M8 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blas95_lp64.a 
# M9 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_cdft_core.so
# M10 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_core.so 
# M11 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_gf_ilp64.so
# M12 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_gf_lp64.so 
# M13 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_gnu_thread.so 
# M14 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_intel_ilp64.so 
# M15 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_intel_thread.so
# M16 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_lapack95_ilp64.a 
# M17 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_lapack95_lp64.a 
# M18 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_pgi_thread.so
# M19 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_scalapack_ilp64.so 
# M20 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_scalapack_lp64.so 
# M21 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_sequential.so
# M22 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_tbb_thread.so
# M23 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_avx.so 
# M24 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_avx2.so 
# M25 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_avx512_mic.so
# M26 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_avx512.so 
# M27 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_def.so 
# M28 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_mc.so
# M29 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_mc3.so 
# M30 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_rt.so 
# M31 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_avx.so 
# M32 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_avx2.so
# M33 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_avx512_mic.so 
# M34 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_avx512.so 
# M35 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_cmpt.so 
# M36 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_def.so 
# M37 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_mc.so 
# M38 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_mc2.so 
# M39 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_vml_mc3.so
# M40 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_pgi_thread.a
# M41 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_scalapack_ilp64.a 
# M42 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_sequential.a
# M43 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_tbb_thread.a
# M44 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_scalapack_lp64.a 
# M45 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_cdft_core.a
# M46 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_core.a 
# M47 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_gf_ilp64.a
# M48 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_gf_lp64.a 
# M49 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_gnu_thread.a 
# M50 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_intel_ilp64.a 
# M51 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_intel_thread.a
# M52 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_intelmpi_ilp64.a
# M53 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_openmpi_ilp64.a
# M54 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_openmpi_lp64.a 
# M55 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_sgimpt_ilp64.a 
# M56 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_sgimpt_lp64.a  
# M57 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_blacs_intelmpi_lp64.so
# M58 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.a
# M59 = /home/stark/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.so
# MKL = $(M46) $(M1) $(M2) $(M3) $(M4) $(M5) $(M6) $(M7) $(M8) $(M9) $(M10) $(M11) $(M12) $(M13) $(M14) $(M15) $(M16) $(M17) $(M18) $(M19) $(M20) $(M21) $(M22) \
# $(M23) $(M24) $(M25) $(M26) $(M27) $(M28) $(M29) $(M30) $(M31) $(M32) $(M33) $(M34) $(M35) $(M36) $(M37) $(M38) $(M39) $(M40) $(M41) $(M42) $(M43) $(M44) \
# $(M45) $(M47) $(M48) $(M49) $(M50) $(M51) $(M52) $(M53) $(M54) $(M55) $(M56) 
# INC = MKL/include 
# PLIB = MKL/intel64_lin
#MKL = -Wl,--start-group $(PLIB)/libmkl_blas95_ilp64.a $(PLIB)/libmkl_lapack95_ilp64.a $(PLIB)/libmkl_intel_ilp64.a $(PLIB)/libmkl_sequential.a $(PLIB)/libmkl_core.a -Wl,--end-group
# MKL = -Wl,--start-group $(PLIB)/libmkl_blas95_ilp64.a $(PLIB)/libmkl_lapack95_ilp64.a $(PLIB)/libmkl_intel_ilp64.a $(PLIB)/libmkl_sequential.a $(PLIB)/libmkl_core.a -Wl,--end-group -ldl
#MKL = -lmkl_blas95_ilp64 -lmkl_lapack95_ilp64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core 


COPT = -O3 -Wall 

default: main

clean: 
	rm *.o 
	rm main

main: main.c prob.o time.o  umfpack.o residue.o plot.o GaussSeidel.o MultiGrid.o SparseTensor.o ConjugateGradient.o 
	$(CC) $(COPT) -g $^ -o $@ $(LIB) -L$(PLIB) -I$(INC) $(MKL) -lpthread -lblas -llapack -lm

umfpack.o: umfpack.c
	$(CC) $(COPT) -g -c $< -o $@ -ISuiteSparse/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include

%.o: %.c
	$(CC) $(COPT) -g -c $< -o $@


