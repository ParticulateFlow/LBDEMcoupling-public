# LBDEMcoupling

* [About](#about)
* [Compatibility](#compatibility)
* [Installation](#installation)
* [Setting Up a Simulation](#setting_up)
* [Implicit Assumptions, Known Issues](#assumptions)
* [Gallery](#gallery)
* [Documented Compiler Switches](#compilerswitches)
* [References](#references)
* [License and Copyright](#license)

<a name="about"></a>
## About

LBDEMcoupling is a coupling between the Lattice Boltzmann (LB) library
Palabos (http://www.palabos.org) and the Discrete Element Method code
LIGGGHTS® (http://www.ligggghts.com). It implements the model of Noble
and Torczinsky [[1]](#ref1) for resolved coupling between particles
and a fluid phase.

<a name="compatibility"></a>
## Compatibility

Currently, LBDEMcoupling is compatible with the following versions:
* LIGGGHTS 3.1
* Palabos v1.5r1

We are working on establishing compatibility with the latest version of LIGGGHTS.

<a name="installation"></a>
## Installation

Just clone the git repo, and create the following three shell
variables in your .bashrc:

          export PALABOS_ROOT=path/to/palabos/
          export LIGGGHTS_ROOT=path/to/liggghts/
          export LBDEM_ROOT=path/to/lbdem/

You also need to move the files

         fix_lb_coupling_onetoone.cpp
         fix_lb_coupling_onetoone.h

to your LIGGGHTS installation directory and recompile LIGGGHTS from
scratch with

         make clean-all; make fedora



<a name="setting_up"></a>
## Setting Up a Simulation

It is assumed that you have basic knowledge in both Palabos and
LIGGGHTS. If you are not familiar with one or both of the codes, it is
strongly recommended that you work through a few tutorial and example
cases.

To compile a case, it is necessary to build LIGGGHTS as a
library. This is explained in the LIGGGHTS documentation. Compiling a
case works just like in Palabos. In the Makefile, you will find the
line

          linkFlags    = -Wl,--whole-archive -llmp_fedora -Wl,--no-whole-archive

In this line, you need to replace the -llmp_fedora with the name of
your library (usually called lmp_whatever).

<a name="gallery"></a>
## Gallery

### Settling spheres

<img src="doc/img/settling.png" alt="10000 settling spheres">

### Square Channel with Particles

<img src="doc/img/showcaseRectChannel.png">


<a name="assumptions"></a>
## Implicit Assumptions, Known Issues

### Sphere size

The code implicitly assumes that all your particles are larger than
four grid spacings and smaller than half the smallest extent of any
partition.

### Multiple Definition Errors

If you come across a multiple definition error during compilation,
please consult ![this
thread](http://www.palabos.org/forum/read.php?11,6746,7581)
in the Palabos forum.

<a name="compilerswitches"></a>
## Documented Compiler Switches

LBDEM_USE_MULTISPHERE switches on support for the (non-public)
multisphere model of LIGGGHTS

<a name="references"></a>
## References

<a name="ref1">[1]</a> Noble, D. R., & Torczynski, J. R. (1998). A
lattice-Boltzmann method for partially saturated computational
cells. *International Journal of Modern Physics C*, 9(08), 1189-1201.

<a name="license"></a>
## License and Copyright

(c) Johannes Kepler University Linz, Austria

released under the GPLv3

main author: Philippe Seil (philippe.seil@jku.at)

LIGGGHTS® is a registered trade mark of DCS Computing GmbH, the
producer of the LIGGGHTS® software
