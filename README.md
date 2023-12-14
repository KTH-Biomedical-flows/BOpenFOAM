# BOpenFOAM
*****

OpenFOAM 10 with an added blood-related solver (vfDependentViscosityFoam). Largely based on the work by [Niclas Berg](https://github.com/niclasberg/BOpenFOAM). The solver considers blood as a single-phase fluid with viscosity that depends on shear rate and haematocrit (different models available). Haematocrit is computed with a convection-diffusion equation (diffusion can be Fickian or Leighton-Acrivos). 
The solver can also be coupled with a LPT function object that computes platelet trajectories and activation both in runtime and as a postprocessing step. An example can be found under OpenFOAM-10/tutorials/multiphase/vfDependentViscosityFoam.

## Build steps

Set appropriate environment variables

```
source OpenFOAM-10/etc/bashrc
```

Compile needed third party code (mainly scotch)

```
cd ThirdParty-10
./Allwmake -j
```

Compile OpenFOAM 

```
cd OpenFOAM-10
./Allwmake -j
```

The ```-j``` flag uses all available cores/hyperthreads (a custom number can also be prescribed).
