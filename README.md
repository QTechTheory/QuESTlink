<img align="left" src="Doc/quest_logo.png" alt="The QuEST logo">
<img align="left" src="Doc/mma_logo.png" alt="The Mathematica logo">

# [QuESTlink](https://questlink.qtechtheory.org/)
A [Mathematica](https://www.wolfram.com/mathematica/) package for remote multithreaded and GPU emulation of quantum computers, using [QuEST](https://quest.qtechtheory.org/).

> QuESTlink is currently considered as an early-release, and is under active development. 
> The offered API may contain bugs and is liable to change.

QuESTlink can be immediately deployed in Mathematica without installation. Simply run
```Mathematica
Import["https://qtechtheory.org/QuESTlink.m"]

CreateDownloadedQuESTEnv["MacOs"]
```

> Pre-prepared Windows and Linux builds of QuESTlink are coming soon.
> In the meantime, read the [Doc](Doc/README.md) to compile QuESTlink from source on these platforms, and/or enable multithreading and GPU-acceleration.

To learn how to use QuESTlink, see our [whitepaper](https://arxiv.org/abs/1912.07904), 
or some demos at [questlink.qtechtheory.org](https://questlink.qtechtheory.org).

If QuESTlink helps you in your research, feel free to cite our paper:
```
Tyson Jones, Simon C Benjamin, 
"QuESTlink - Mathematica embiggened by a hardware-optimised quantum emulator"
arXiv:1912.07904 (2019)
```
