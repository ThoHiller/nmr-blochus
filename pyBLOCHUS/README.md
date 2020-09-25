## pyBLOCHUS

Python implementation of the **BLOCHUS** core functionality without any graphical user interface.

### Table of Contents
1. [About](#about)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Usage](#usage)

- - -
<a name="about"></a>
### About

**pyBLOCHUS** is a set of python modules, that allow some basic simulations of (S)NMR spin dynamics based on the Bloch equations. The Bloch equations are solved in the laboratory frame of reference with the Scipy solve-ip solver. Because it was developed within the scope of a near surface SNMR project, its main features are the simulation of (1) pre-polarization switch-off ramps and (2) excitation pulses. For graphics output the [matplotlib](https://matplotlib.org/) library is used.

- - -
<a name="requirements"></a>
### Requirements

In order to work properly you need to meet the following requirements:

1. A working [Python](https://www.python.org/) installation (I recommend [Anaconda](https://www.anaconda.com/products/individual) - as this is what I use under Windows)
2. The [NumPy](https://numpy.org/), [SciPy](https://www.scipy.org/) and [matplotlib](https://matplotlib.org/) libraries are <span style="color:red">required</span>
2. The [Numba](https://numba.pydata.org/) JIT compiler is optional but speeds up the computations considerably.

#### Operating System

I tested it successfully under Windows 7 (64bit) and 10 (64bit) with [Anaconda](https://www.anaconda.com/products/individual) having installed [Python](https://www.python.org/) 3.7.7. and the [Spyder](https://www.spyder-ide.org/) IDE 4.1.3. .

**NOTE:** So far I did not test anything on Linux or a Mac. If you get it to work on either of the two systems (which it basically should I guess) please let me know.

- - -
<a name="installation"></a>
### Installation

1. Makes sure that all required packages are installed and work properly
2. Download or `clone` the **BLOCHUS** repository and put the folder somewhere on your hard drive.
3. As I recommend uying an IDE like [Spyder](https://www.spyder-ide.org/), make sure that the **pyBLOCHUS** folder is on your PYTHONPATH (e.g. in [Spyder](https://www.spyder-ide.org/):  Tools -> PYTHONPATH manager -> Add path)
4. Sometimes a [Spyder](https://www.spyder-ide.org/) restart is necessary tu update/sync the PYTHONPATH

- - -
<a name="usage"></a>
### Usage

Have a look in the [example](pyBLOCHUS/examples) scripts  for the usage of the different modules. 
