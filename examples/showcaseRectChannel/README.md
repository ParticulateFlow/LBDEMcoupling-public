# ShowcaseRectChannel

In this case, flow through a square channel with particles is
simualated. 

<img src="../../doc/img/showcaseRectChannel.png">

The channel crossection is 0.2x0.2m. Viscosity is hard-wired to 1e-3
m^2/s, density is set to 1000 kg/m^3. All other parameters can be set
via the command line. The flow is driven by a periodic pressure
boundary condition by Zhang and Kwok ![[1]](#ref1)

## Running the case

The command line for the case reads

```
./showcaseRectChannel N deltaP v_frac d_oart uMax lChannel time outDir
```

with the parameters
* N being the number of grid points along one side of the crossection
* deltaP being the pressure difference between in- and outlet in Pa
* v_frac being the particle volume fraction
* d_part being the particle diameter
* uMax being the maximum LB velocity, proportional to the Mach number
* lChannel being the length of the channel as a multiple of the width
* time being the runtime as multiple of the Stokes time
* outDir being the directory for your output files. The out dir must
contain subdirectories outDir/post and outDir/tmp for LIGGGHTS and
Palabos data respectively.

## Results

Depending on solid fraction and Reynolds number, pattern formation
(clusters, rings) should occur, similar to the Segre-Silberberg effect
![[2]](#ref2). However, this can take O(10000) times
the Stokes time, so please be patient and let the simulation run.


<a name="ref1">[1]</a> Zhang, J., & Kwok, D. Y. (2006). Pressure boundary
condition of the lattice Boltzmann method for fully developed periodic
flows. *Physical review E*, 73(4), 047702.

<a name="ref2">[2]</a> Segre, G., & Silberberg, A. (1962). Behaviour
of macroscopic rigid spheres in Poiseuille flow Part 2. Experimental
results and interpretation. *Journal of fluid mechanics*, 14(01),
136-157.