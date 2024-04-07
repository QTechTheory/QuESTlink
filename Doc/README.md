
# Doc

If you wish to run QuESTlink in multithreaded or GPU-accelerated modes, or as a server, you will need to download and manually compile it. Below are instructions for getting QuESTlink up and running, when it doesn't quite run out-of-the-box.

  * [Tools](#tools)
  * [Download](#download)
  * [Compile](#compile)
  * [Launch offline](#launch-offline)
  * [Launch as a server](#launch-as-a-server)

______________


## Tools 

Don't have a compiler handy? We recommend the [comprehensive guide](https://quest.qtechtheory.org/download/) for the parent QuEST project to obtain the tools needed below.


## Download

The best way to download the QuESTlink source code is through [`git`](https://git-scm.com/) at the command-line.

```bash
git clone --recurse-submodules https://github.com/QTechTheory/QuESTlink.git
cd QuESTlink
```
This will also download the [QuEST](https://github.com/QuEST-Kit/QuEST) submodule source-code.
You are thereafter in the root directory of QuESTlink and ready to compile.

> If you wish to compile a specific *Github branch* like `develop`, simply additionally run
> ```bash
> git checkout develop
> ```
> You can also perform this within the `/QuEST` subdirectory to change the branch of the QuEST submodule and access in-development features.

______________________________

## Compile

Compiling is trivial with [GNUMake](https://www.gnu.org/software/make/) and the provided [makefile](../makefile), and a C++ compiler (which supports `C++11`).

> See [here](WINDOWS.md) for a comprehensive guide to compiling QuESTlink on **Windows**, including how to obtain the necessary compilers.

> Note **Linux** users should first run `sudo apt-get install uuid-dev` before compiling. If this fails (and/or compiling says `cannot find -luuid`), also try to install `uuid`, `uuid-devel` and `libuuid-devel`. On a SLURM cluster, this can often be resolved with `module load util-linux`.

Within the root directory, edit the [makefile](../makefile) and set:

- `OS` to your operating system (`LINUX`, `MACOS` or `WINDOWS`)

- `COMPILER` to your C++11 compiler command. 
  > If in doubt, leave this as `g++`
- `COMPILER_TYPE` to the type of your compiler (`GNU`, `INTEL`, `CLANG` or `MSVC`).
  > If in doubt, run `g++ --version` in terminal for a clue. Otherwise, `COMPILER_TYPE` will likely match your `OS`, as: `LINUX` & `GNU`, `MACOS` & `CLANG`, `WINDOWS` & `MSVC`.
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

## Launch offline 

QuESTlink can be launched without an internet connection, using a copy of this repository. In Mathematica, simply run 

```Mathematica 
Import["path/to/QuESTlink/Link/QuESlink.m"];
CreateLocalQuESTEnv["path/to/quest_link"];
```

where `quest_link` has been compiled as above, or previously obtained using 
```Mathematica 
CreateDownloadedQuESTEnv[];
```

## Launch as a server

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
