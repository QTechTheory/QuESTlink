## remove all current object files
rm *.o
rm quest_link

## compile QUEST backend
echo "Compiling QuEST backend..."
gcc-8 -x c -O2 -std=c99 -mavx -Wall -fopenmp -DQuEST_PREC=2  -c QuEST/QuEST.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -fopenmp -DQuEST_PREC=2  -c QuEST/QuEST_validation.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -fopenmp -DQuEST_PREC=2  -c QuEST/QuEST_common.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -fopenmp -DQuEST_PREC=2  -c QuEST/QuEST_qasm.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -fopenmp -DQuEST_PREC=2  -c QuEST/mt19937ar.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -fopenmp -DQuEST_PREC=2  -c QuEST/CPU/QuEST_cpu.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -fopenmp -DQuEST_PREC=2  -c QuEST/CPU/QuEST_cpu_local.c

## compile MMA functions
echo "Compiling Mathematica wrappers..."
gcc-8 -x c -O2 -std=c99 -mavx -Wall -fopenmp -arch x86_64 -DQuEST_PREC=2  -IQuEST/CPU -IQuEST -c quest_link.c

## generate template sources
echo "Generating template sources..."
./WSTPlibs/wsprep quest_templates.tm -o quest_templates.tm.c
echo "Compiling template sources..."
gcc-8 -fopenmp -arch x86_64 -c -o quest_templates.o quest_templates.tm.c

## link
echo "Linking QuEST, wrappers and templates..."
gcc-8 -O2 -std=c99 -mavx -Wall -fopenmp -arch x86_64 -DQuEST_PREC=2  -IQuEST/CPU -IQuEST -o quest_link quest_link.o quest_templates.o QuEST.o QuEST_validation.o QuEST_common.o QuEST_qasm.o mt19937ar.o QuEST_cpu.o QuEST_cpu_local.o -lm -lc++ WSTPlibs/libWSTPi4.36.a -framework Foundation

# cleanup
echo "Cleaning up..."
##rm quest_templates.tm.c
rm *.o

echo "Done!"
