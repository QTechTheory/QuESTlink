
# Windows

Below are instructions for compiling QuESTlink from source on Windows, using only the command line. These instructions will install the Microsoft Visual Studio Compiler (MSVC), the Chocolatey package manager, and a Windows-compatible GNUMake.

- [1 - Build Tools for VS](#1---build-tools-for-vs)
- [2 - Choco](#2---choco)
- [3 - GNUMake](#3---gnumake)
- [4 - Compile](#4---compile)
- [5 - Run](#5---run)

### 1 - Build Tools for VS

Download **Build Tools for Visual Studio**  from [this page](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019), and install it. 

This will provide both the **Developer Powershell for VS** and the **Developer Command Prompt for VS**.

Open **Developer Command Prompt** and check that commands

```bash 
cl
```
and 
```bash 
link -help
```
are correctly recognised.

### 2 - Choco 

Open **Developer Powershell**  as an administrator and enter 
```bash 
Set-ExecutionPolicy AllSigned
``` 
then
```bash 
Set-ExecutionPolicy Bypass -Scope Process -Force; [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072; iex ((New-Object System.Net.WebClient).DownloadString('https://chocolatey.org/install.ps1'))
```

> If you are unable to paste this command via ctrl-v, click the icon top-left of the powershell window and select paste

This downloads [Chocolatey](https://chocolatey.org/), a Windows package manager.

### 3 - GNUMake 

Still in **Developer Powershell**, enter 
```bash 
choco install make
```
This creates the `make` command in the **Developer Command Prompt**.

### 4 - Compile 

We're now ready to compile. Open **Developer Command Prompt** and make sure 
```bash 
make --version
```
is a recognised command; you may have to close and re-open the command prompt first. 

Next, open the [`makefile`](../makefile) in any editor, and set:
- `OS = WINDOWS`
- `COMPILER = cl`
- `COMPILER_TYPE = MSVC`
- `WINDOWS_ARCH = 64` if using 64-bit Windows, else `32` (for x86)

> If using 64-bit Windows, you must use a *64-bit* **Developer Command Prompt** for the following instructions. E.g. `VS2019 x64 Native Tools Command Prompt`. You can find this prompt in the same directory as the developer command prompt. 

To compile for GPU mode, install the [CUDA toolkit](https://developer.nvidia.com/cuda-downloads) (to get the `nvcc` command), and additionally set
- `GPUACCELERATED = 1`
- `GPU_COMPUTE_CAPABILITY = ` to the value corresponding to your GPU (look up [here](https://developer.nvidia.com/cuda-gpus)) with no decimal-point (e.g. `7.1` becomes `71`).

Then, in the **Developer Command Prompt**, navigate to the root QuESTlink directory (where [`makefile`](../makefile) is located) and run 
```bash 
make
```
If successful, the `quest_link.exe` executable will be created, along with several leftover `.o` files which can be safely removed with 
```bash 
make tidy 
```
To recompile after changing a setting in the makefile, run
```bash
make clean
make
```

> Note that trying to run `quest_link.exe` directly at this stage will report a DLL error. 

### 5 - Run 

#### 5.1 - Locally

You can now open Mathematica and run 

```Mathematica 
SetDirectory["path/to/QuESTlink/"]

Import["Link/QuESTlink.m"]
CreateLocalQuESTEnv["quest_link.exe"]
```
and use all facilities of QuEST.



#### 5.2 - Remotely

The created `quest_link.exe` can be used as a server, accessed by Mathematica on another machine. To do this, after compiling, you must copy one of the following DLLs (depending on whether `WINDOWS_ARCH` is 32 or 64 bit) to the same location as `quest_link.exe`.

```bash 
copy WSTP\Windows\wstp32i4.dll .
copy WSTP\Windows\wstp64i4.dll .
```
Running `quest_link.exe` directly now will create a network prompt, which can be ignored/closed.

To launch the server from the **Developer Command Prompt**, run
```bash
quest_link.exe -linkcreate -linkprotocol TCPIP -linkname <PORT1>@<IP>,<PORT2>@<IP>
```
substituting <PORT1> and <PORT2> with two open and available ports, and <IP> with the server IP or domain name.

Then in your Mathematica kernel on another machine, connect to it via

```Mathematica
CreateRemoteQuESTEnv[<IP>, <PORT1>, <PORT2>];
```
