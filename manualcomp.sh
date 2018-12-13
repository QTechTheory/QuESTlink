
## compile QUEST backend
echo "Compiling QuEST backend..."
gcc-8 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/QuEST.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/QuEST_validation.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/QuEST_common.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/QuEST_qasm.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/mt19937ar.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/CPU/QuEST_cpu.c
gcc-8 -x c -O2 -std=c99 -mavx -Wall -DQuEST_PREC=2  -c QuEST/CPU/QuEST_cpu_local.c

## compile MMA functions
echo "Compiling Mathematica wrappers..."
gcc-8 -x c -O2 -std=c99 -mavx -Wall -arch x86_64 -DQuEST_PREC=2  -IQuEST/CPU -IQuEST -c quest_wrapper.c

## generate template sources
echo "Generating template sources..."
./WSTPlibs/wsprep mytest.tm -o mytest.tm.c
echo "Compiling template sources..."
gcc-8 -arch x86_64 -c -o mytest.o mytest.tm.c

## link
echo "Linking QuEST, wrappers and templates..."
gcc-8 -O2 -std=c99 -mavx -Wall -arch x86_64 -DQuEST_PREC=2  -IQuEST/CPU -IQuEST -o quest_link mytest.o QuEST.o QuEST_validation.o QuEST_common.o QuEST_qasm.o mt19937ar.o QuEST_cpu.o QuEST_cpu_local.o quest_wrapper.o -lm -lc++ WSTPlibs/libWSTPi4.36.a -framework Foundation

# cleanup
echo "Cleaning up..."
rm mytest.tm.c
rm *.o

echo "Done!"