#----------------------------------------------
# makefile BGK3D_GPU
# Giorgio Amati, CASPUR 2000-2012
# Giorgio Amati, CINECA 2013-202?

# targets:
#
#		<none> default ---> GPU code with do concurrent (nvfortran)
#		offload ---> GPU code with openmp offload
#		openacc ---> GPU code with openACC
#		serial  ---> Serial CPU code
#		mulitcore -> Multicore CPU code 
# 
#----------------------------------------------
DATE = $(shell date +%d%m%y)
FINE = "******** End Compilation step  *******"
CHECK1 = NVIDIA
CHECK2 = SINGLE
#----------------------------------------------
# Compilers & Flags
#----------------------------------------------
# flags
# default (GNU)
MULTI =   -fopenmp -DMULTICORE -ftree-parallelize-loops=4
OPENACC = -fopenacc -DNOMANAGED
OFFLOAD = -fopenmp  -foffload=nvptx-none -DOFFLOAD
DC = not supported
SFC = gfortran
SCC = gcc
FOPT =  -Ofast
COPT =  -O3
OPT =   -Ofast
DEBUG = -Wuninitialized -DDEBUG_1 -Wunused
#
CHECK1=GNU
CHECK2=SINGLE
CHECK3=ORIGINAL
CHECK4=LDC
#
ifdef NVIDIA
 OFFLOAD = -mp=gpu -DOFFLOAD  -Minfo=accel 
 OPENACC = -Minfo=acc -O2 -acc -DOPENACC -ta=tesla:cc80,managed
 #OPENACC = -Minfo=acc -O2 -acc -DOPENACC -ta=tesla:cc70,managed
 MULTI = -stdpar=multicore -Minfo=all -DMULTICORE
 DC  = -stdpar=gpu -DDC
 SCC = nvcc 
 SFC = scorep-nvfortran  -DPGI
 SFC = nvfortran -DPGI
 FOPT = -O2 -Mnodepchk -Mcontiguous
 COPT = -O2
 OPT  = -O2 -Mnodepchk -Mcontiguous
 CHECK1=NVIDIA
endif
#
ifdef AMD
 MULTI = -DMULTICORE
 OFFLOAD = -fopenmp --offload-arch=gfx90a -DOFFLOAD
 OPENACC = not supported
 DC = not supported
 SFC = flang -DPGI
 SCC = clang -DPGI
 FOPT =  -O3
 COPT =  -O3
 OPT =   -O3
 CHECK1=AMD
endif
#
ifdef CRAY
 MULTI =
 OFFLOAD = -DPGI -h omp -DOFFLOAD -DCRAY
 OPENACC = -DPGI -h acc -DNOMANAGED -DOPENACC
 DC = not supported
 SCC = cc
 SFC = ftn
 FOPT = -O3
 COPT =
 OPT = -O3
 CHECK1=CRAY
endif
#
ifdef GNU
# do nothing (default)	
endif
#
ifdef INTEL
 MULTI =  -qopenmp -DMULTICORE
 OPENACC = not supported
 OFFLOAD = 
 DC =  -DDC
 SFC = ifx
 SCC = icx
 FOPT = -O3 -xCORE-AVX512 -mtune=skylake-avx512 -assume contiguous_pointer
 COPT = -O3 
 OPT =  -O3
 OPENMP = -qopenmp
 INFO = -warn all
 CHECK1=INTEL
endif
#
#
# to remove
#ifdef PWR
# SCC = xlc_r  
# SFC = xlf_r 
# OFFLOAD = -qsmp=omp -qoffload -DOFFLOAD
# FOPT = -O3 -qassert=contig
# COPT = -O3 -qassert=contig
# OPT = -O3 -qassert=contig
# CHECK1=PWR
#endif
#
#----------------------------------------------
# precision used
#----------------------------------------------
#
# default
ifdef SINGLE
 PREC =
 CHECK2=SINGLE
endif
#
ifdef MIXED1
 PREC =-DHALF_P -DMIXEDPRECISION 
 CHECK2=MIXED1
endif
#
ifdef MIXED2
 PREC =-DMIXEDPRECISION 
 CHECK2=MIXED2
endif
#
ifdef DOUBLE
 PREC =-DDOUBLE_P
 CHECK2=DOUBLE
endif
#
#----------------------------------------------
# implementations 
#----------------------------------------------
#
# default, ORIGINAL
ifdef ORIG
 IMPL = 
 CHECK3=ORIG
endif
#
ifdef FUSED
 IMPL = -DFUSED
 CHECK3=FUSED
endif
#
#
#----------------------------------------------
# test case
#----------------------------------------------
#
# default, Lid Driven Cavity
ifdef LDC
 PREPROC = -DLDC
 CHECK4=LDC
endif
#
# Poiseulle Flow
ifdef POF
 PREPROC = -DCHANNEL -DPOF
 CHECK4=POF
endif
#
# Taylor Green Vortex
ifdef TGV
 PREPROC = -DPERIODIC -DTGV
 CHECK4=TGV
endif
#
# Von Karman Streets  (SPHERE)
ifdef VKS
 PREPROC =-DINFLOW -DOBSTACLE -DDRAG -DVKS
 CHECK4=VKS
endif
#
# Kida Vortex
ifdef KVX
 PREPROC =-DPERIODIC -DKVX
 CHECK4=KVX
endif
#
# Homogenous Isotropic Turbulence  (not forced)
ifdef HIT
 PREPROC =-DPERIODIC -DTGV3D -DHIT
 PREPROC =-DPERIODIC -DTGV3D -DVTK3D -DHIT
 CHECK4=HIT
endif
#
#----------------------------------------------
# module & subroutines
#----------------------------------------------
MOD0 = get_mem_mod.o \

MOD1 = precision.o \
      storage.o \
      mod_tools.o \

# 0th level subroutine
OBJ0 = bgk3d.o \

# 1st level subroutine
OBJ1 = setup.o \
       initialize.o \
       diagnostic.o \
       boundaries.o \
       propagation.o \
       save.o\
       finalize.o \

# 2nd level subroutine
OBJ2 = check_case.o \
       input.o \
       alloca.o \
       bcond.driven.o \
       bcond.periodic.o \
       bcond.channel.o \
       bcond.inflow.o \
       bcond.obs.o \
       build.bcond.o \
       build.obs.o \
       restore_raw.o \
       init.o\
       hencol.o \
       outdat.o \
       vtk_obstacle.o \
       vtk_3d.binary.o \
       vtk_xy.o \
       vtk_xy.binary.o \
       vtk_xz.binary.o \
       vtk_yz.binary.o \
       vtk_om.xy.binary.o \
       varm.o \
       prof_k.o \
       prof_j.o \
       prof_i.o \
       diagno.o \
       dissipation.o \
       dissipation2.o \
       probe.o \
       movef.o\
       save_raw.o\
       save_vel.o\
       profile.o\
       git_info.o \
       make_info.o \
       draglift.o \
       dealloca.o

# 3rd level subroutine
OBJ3 = get_mem.o \
       collision.o \
       col_MC.o \
       col.o 

OBJ = $(OBJ1) $(OBJ0) $(OBJ2) $(OBJ3)

#----------------------------------------------
# targets
#----------------------------------------------
#
#
serial: all
#
#----------------------------------------------
all: VER := serial
all: FC := $(SFC) 
all: CC := $(SCC)
all: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(LIB) $(FIX) -o ../RUN/bgk3d.$(VER).x 
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4) 
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
doconcurrent: VER := doconcurrent
doconcurrent: FOPT += $(DC)
doconcurrent: OPT  += $(DC)
doconcurrent: FC := $(SFC) 
doconcurrent: CC := $(SCC)
doconcurrent: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(LIB) $(FIX) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
offload: FOPT += $(OFFLOAD)
offload: OPT  += $(OFFLOAD)
offload: VER := offload
offload: FC := $(SFC)
offload: CC := $(SCC)
offload:$(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(FIX) $(LIB) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4) 
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
openacc: FOPT += $(OPENACC)
openacc: OPT  += $(OPENACC)
openacc: VER := openacc
openacc: FC := $(SFC)
openacc: CC := $(SCC)
openacc:$(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(FIX) $(LIB) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4) 
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
multicore: VER:= multicore
multicore: FC := $(SFC) $(MULTI)
multicore: CC := $(SCC)
multicore: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(FIX) $(LIB) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
# rules 
#----------------------------------------------

%.o %.mod: $(INC) %.f90
	$(FC) $(FOPT) $(PREC) $(PREPROC) $(IMPL) $(FIX) -c $<

%.o: $(INC) %.F90
	$(FC) $(FOPT) $(PREC) $(PREPROC) $(IMPL) $(FIX) -c $<

%.o: $(INC) %.f
	$(FC) $(FOPT) $(PREC) $(PREPROC) $(IMPL) $(FIX) -c $<

%.o: $(INC) %.c
	$(CC) $(COPT)                          -c $<

make_info.F90:
	echo "    subroutine make_info"                             >  make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(CHECK1) ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(CHECK2) ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(CHECK3) ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(CHECK4) ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FIX)    ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FC)     ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FOPT)   ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(PREPROC)` \"         >> make_info.F90
	echo "    end subroutine make_info"                         >> make_info.F90
	
git_info.F90:
	echo "    subroutine git_info"                              >> git_info.F90
	echo "    write( 6,*) \""INFO: `git log  | head  -n 1` \"   >> git_info.F90
	echo "    write(16,*) \""INFO: `git log  | head  -n 1` \"   >> git_info.F90
	echo "    end subroutine git_info"                          >> git_info.F90
#
#----------------------------------------------
# utilities
#----------------------------------------------
clean:
	rm -f *.o;
	rm -f *.mod;
	rm -f core;
#----------------------------------------------
info:  $(INC)
	@echo "Modules          =  "$(MOD0) $(MOD1);
	@echo "Objects          =  "$(OBJ);
	@echo "Compiler         =  "$(SFC);
	@echo "Compiler  flags  =  "$(OPT);
	@echo "Preproc   flags  =  "$(PREPROC);
	@echo "Fix       flags  =  "$(FIX);
	@echo "Precision flags  =  "$(CHECK2);
	@echo "Implemet. flags  =  "$(CHECK3);
	@echo "Version          =  "$(CHECK4);

#----------------------------------------------
test:
	make clean;
	make serial;
	make clean;
	make multicore;
	make clean;
	make doconcurrent;
	make clean;
	make openacc;
	make clean;
	make offload;
#----------------------------------------------
tar:
	make clean
	rm -rf ../RUN/*.vtk
	rm -rf ../RUN/*.bin
	rm -rf ../RUN/*.h5
#	cd ../TEST/
#	./clean.x
#	cd -
	tar -cvf ../../archive.bgk3d.$(DATE).tar ../../BGK_3D/*
#----------------------------------------------
help:
	@echo "different target available"
	@echo "do concurrent (default, GPU) "
	@echo "serial (CPU)"
	@echo "multicore (CPU)"
	@echo "offload (GPU)"
	@echo "openacc (GPU)"
#----------------------------------------------

