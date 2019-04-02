
rm quest_link

## compile QUEST backend
echo "Compiling QuEST backend..."
clang-3.7 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/QuEST.c
clang-3.7 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/QuEST_validation.c
clang-3.7 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/QuEST_common.c
clang-3.7 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/QuEST_qasm.c
clang-3.7 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/mt19937ar.c

nvcc -dc -O2 -arch=compute_61 -code=sm_61 -DQuEST_PREC=2 -ccbin clang-3.7 QuEST/GPU/QuEST_gpu.cu
##clang-3.9 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/CPU/QuEST_cpu.c
##clang-3.9 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/CPU/QuEST_cpu_local.c

## compile MMA functions
echo "Compiling Mathematica wrappers..."
clang-3.7 -x c -O2 -std=c99 -mavx -Wall -arch x86_64 -DQuEST_PREC=2  -IQuEST/CPU -IQuEST -c quest_link.c

## generate template sources
echo "Generating template sources..."
./WSTPlibs/wsprep quest_templates.tm -o quest_templates.tm.c
echo "Compiling template sources..."
clang-3.7 -arch x86_64 -c -o quest_templates.o quest_templates.tm.c

## link
echo "Linking QuEST, wrappers and templates..."
nvcc -O2 -arch=compute_61 -code=sm_61 -DQuEST_PREC=2 -ccbin clang-3.7 -IQuEST/GPU -IQuEST -o quest_link quest_link.o quest_templates.o QuEST.o QuEST_validation.o QuEST_common.o QuEST_qasm.o mt19937ar.o QuEST_gpu.o -lm -lc++ WSTPlibs/libWSTPi4.36.a ## -framework Foundation
## gcc-8 -O2 -std=c99 -mavx -Wall -arch x86_64 -DQuEST_PREC=2  -IQuEST/CPU -IQuEST -o quest_link quest_link.o quest_templates.o QuEST.o QuEST_validation.o QuEST_common.o QuEST_qasm.o mt19937ar.o QuEST_cpu.o QuEST_cpu_local.o -lm -lc++ WSTPlibs/libWSTPi4.36.a -framework Foundation



# cleanup
echo "Cleaning up..."
rm quest_templates.tm.c
rm *.o

echo "Done!"
