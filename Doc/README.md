
# Doc

Our [whitepaper](https://arxiv.org/abs/1912.07904) gives a thorough overview of the QuEST syntax and functionality. This, and our [demos](../Demos), are the best ways to learn QuESTlink. Learn more at [quest.qtechtheory.org](https://questlink.qtechtheory.org/). 

Below are instructions for getting QuESTlink up and running, when it doesn't quite run out-of-the-box.

  * [Running](#running)
  * [Compiling](#compiling)
  * [Launching offline](#launching-offline)
  * [Launching as server](#launching-as-server)

______________

## Running

Single-thread QuESTlink can be used immediately from within Mathematica:
```Mathematica 
Import["https://qtechtheory.org/questlink.m"]
CreateDownloadedQuESTEnv[]
```

These commands download the QuESTlink Mathematica package file, and the `quest_link` executable, to your machine.

To see available all available functions and circuit symbols, run
```Mathematica
?QuEST`*
```

______________________________

## Compiling

Multithreaded and GPU-accelerated modes require compiling QuESTlink from source, using [GNUMake](https://www.gnu.org/software/make/) with the provided [makefile](../makefile), and a C++ compiler (which supports C++11).

> See [here](WINDOWS.md) for complete guide on compiling on **Windows**, including how to obtain the necessary compilers.

> Note **Linux** users should first run `sudo apt-get install uuid-dev` before compiling. If this fails (and/or compiling says `cannot find -luuid`, also try to install `uuid`, `uuid-devel` and `libuuid-devel`.

Edit [makefile](../makefile) and set:

- `OS` to your operating system (`LINUX`, `MACOS` or `WINDOWS`)

- `COMPILER` to your C++11 compiler command. 
>If in doubt, leave this as `g++`
- `COMPILER_TYPE` to the type of your compiler (`GNU`, `INTEL`, `CLANG` or `MSVC`).
>If in doubt, run `g++ --version` in terminal for a clue. Otherwise, `COMPILER_TYPE` will likely match your `OS`, as: `LINUX` & `GNU`, `MACOS` & `CLANG`, `WINDOWS` & `MSVC`.
- `MULTITHREADED = 1` to compile in multithreaded mode
> After compiling, and before calling `CreateLocalQuESTEnv[]`, you should set (in terminal)
> ```bash 
> export OMP_NUM_THREADS=<nthreads>
> ```
> replacing `<nthreads>` with the number of CPU cores of your machine, minus a few (sparing them for the Mathematica kernel itself).
>
> **Note** `MULTITHREADED` mode requires an [OpenMP](https://scc.ustc.edu.cn/zlsc/sugon/intel/compiler_f/main_for/optaps/common/optaps_par_openmp_multiple_compilers.htm)-compatible compiler. MacOS users should thus avoid `clang`, and download/use a GNU compiler (e.g. [`gcc@8`](https://formulae.brew.sh/formula/gcc@8)).
- `GPUACCELERATED = 1` to use an NVIDIA GPU.
> This requires an [`nvcc`](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html) compiler compatible with your chosen `COMPILER`. The `nvcc` compiler command can be changed by overwriting `CUDA_COMPILER` in the makefile.
> An `nvcc` compiler can be obtained on Linux with `sudo apt install nvidia-cuda-toolkit`
>
> **Note** you must also set `GPU_COMPUTE_CAPABILITY` in the makefile to the CC corresponding to your GPU. You can look this up [here](https://developer.nvidia.com/cuda-gpus).


With these settings set, QuESTlink is compiled from terminal, in the root directory [`QuESTlink/`](../), via
```bash
make
```

Compiling will create an executable `quest_link`, and some other build files. These other unneeded files can be removed with 
```bash
make tidy 
```
leaving only `quest_link`.

Should you wish to remove all compiled files and start again, run

```bash 
make clean 
```

From within Mathematica, the compiled `quest_link` environment is connected to via 
```Mathematica 
Import[...]
CreateLocalQuESTEnv["path/to/quest_link"];
```

_______________________________

## Launching offline 

QuESTlink can be launched without an internet connection, using a copy of this repository. In Mathematica, simply run 

```Mathematica 
Import["path/to/QuESTlink/Link/QuESlink.m"];
CreateLocalQuESTEnv["path/to/quest_link"];
```

where `quest_link` has been compiled as above, or previously obtained using 
```Mathematica 
CreateDownloadedQuESTEnv[];
```

## Launching as server

To launch `quest_link` as a server, to access remotely from a local kernel, run
```bash
./quest_link -linkcreate -linkprotocol TCPIP -linkname <PORT1>@<IP>,<PORT2>@<IP>
```
substituting `<PORT1>` and `<PORT2>` with two open and available ports, and
`<IP>` with the server IP or domain name.

Then in your local Mathematica kernel, connect to it via 
```Mathematica 
CreateRemoteQuESTEnv[<IP>, <PORT1>, <PORT2>];
```
