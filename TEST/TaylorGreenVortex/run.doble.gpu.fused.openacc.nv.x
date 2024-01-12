#!/bin/tcsh
#
setenv DIR RUN_DOUBLE_GPU_FUSED_OPENACC_NV
setenv EXE bgk3d.openacc.x
#
echo "---------------------------"
echo "starting test TG vortices  "
echo " ---> nvfortran            "
echo " ---> double precision     "
echo " ---> fused                "
echo " ---> openacc              "
echo " ---> " $EXE
echo " ---> " $DIR
echo "---------------------------"
#
rm -r $DIR
mkdir $DIR
cd $DIR

# step 1: compiling
echo "step 1: compiling"
cd ../../../SRC
make clean
make openacc NVIDIA=1 DOUBLE=1 TGV=1 FUSED=1
if ($?) then
   echo "compiling fails..."
   exit 1
else
   cd -
   cp ../../../RUN/$EXE  .
   cp ../bgk.input .
   echo "compiling  ended succesfully..."
endif


# step 2: running test
echo "step 2: running test"
./$EXE >& out.log
if (-e "bgk.perf") then
   echo "run ended succesfully..."
else
   echo "running test fails..."
   exit 1
endif

