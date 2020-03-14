
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
cl --version
```
and 
```bash 
link --version
```
are correctly recognised.

### 2 - Choco 

Open **Developer Powershell** and enter 
```bash 
Set-ExecutionPolicy AllSigned
``` 
then
```bash 
Set-ExecutionPolicy Bypass -Scope Process -Force; [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072; iex ((New-Object System.Net.WebClient).DownloadString('https://chocolatey.org/install.ps1'))
```
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

Then, in the **Developer Command Prompt**, navigate to the root QuESTlink directory (where [`makefile`](../makefile) is located) and run 
```bash 
make
```
If successful, the `quest_link.exe` executable will be created, along with several leftover `.o` files which can be safely removed with 
```bash 
make tidy 
```
> Trying to run `quest_link.exe` directly at this stage will report a DLL error. 

Next, run 
```bash 
copy WSTP\Windws\wstp32i4.dll .
copy WSTP\Windws\wstp64i4.dll .
```
This copies the needed `.dll` files to the same location as `quest_link.exe`. Running `quest_link.exe` directly now will create a network prompt, which can be ignored/closed.

### 5 - Run 

With `wstp32i4.dll` (or `wstp64i4.dll`) in the same location as `quest_link.exe`, you can now open Mathematica and run 

```Mathematica 
SetDirectory["path/to/QuESTlink/"]

Import["Link/QuESTlink.m"]
CreateLocalQuESTEnv["quest_link.exe"]
```
and use all facilities of QuEST.