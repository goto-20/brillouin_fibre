# brillouin_fibre

## Acousto-optics in Optical Waveguides under Mechanical Stress

[brillouin_fibre](https://github.com/goto-20/brillouin_fibre) is a simulation environment for calculating Stimulated Brillouin Scattering (SBS) within a 2-D Optical Waveguide in the presence of external mechanical load. SBS can be described as the process of exciting acoustic waves within a waveguide using two counter-propagating waves.

![SBS](images/BOTDA.png "Stimulated Brillouin Scattering (SBS)")

[brillouin_fibre](https://github.com/goto-20/brillouin_fibre) combines the linear elastic equation, Maxwell's equations and the elastodynamic equations to fully characterize optical waveguides built from any material. The sequence can be described with the following workflow:

![Acousto-optics in Optical Waveguides under Mechanical Stress](images/SBS_WorkFlow.png "Acousto-optics in Optical Waveguides under Mechanical Stress")

The simulation environment is built over the [FEniCS](https://fenicsproject.org/) computing platform for solving partial differential equations (PDEs) with the finite element method (FEM).

## Examples

### Single-Mode Fiber (SMF-28) (Silica)

We make sure that the SBS process for SMF-28 is in agreement with known experimental values.
