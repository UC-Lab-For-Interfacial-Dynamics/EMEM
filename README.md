# Embedded Multiscale Evaporation Model (EMEM)

## Introduction
EMEM is a tuning coefficient-free kinetic evaporation model for the extended meniscus region. EMEM is provided as a set of User Defined Functions (UDFs) to be embedded inside an Ansys Fluent simulation. 

## Theoretical 
The development of this model is motivated by the multiscale nature of the menisci of wetting fluids, where significant evaporation flux in contributed by the transition thin film region and cannot be resolved within macroscale CFD simulations. EMEM uses a differential form of the augmented Young-Laplace equation to construct an ODE describing the evolution of the iterface shape. The ODE is solved numerically using a fixed-step fourth-order Runge-Kutta method at each iteration of the CFD simulation. The Navier-Stokes equations are reduced under the lubrication approximation and solved to obtain the hydrodynamic pressure. One-dimensional heat conduction from the solid wall to the interface is coupled with a latent heat loss condition based on kinetic evaporation. A two-way coupling is established with the CFD simulation: 
    (1) solid wall temperatures are extracted and used as boundary conditions at each CFD iteration
    (2) source terms are returned to the CFD simulation to add the total evaporating mass as a source term at the average-interface temperature.

The code provided here is set up to simulate the extended meniscus of liquid Hydrogen in a 10 mm cylindrical container for validation with the experimental test case.

## References
See associated paper for details:
    (1) MDPI fuels
    (2) Yasin, A.; Bellur, K. Computational Modeling of Evaporation without Tuning Coefficients. Applied Thermal Engineering 2025, 276, 126807. https://doi.org/10.1016/j.applthermaleng.2025.126807.
    (3) Bellur, K.; Médici, E. F.; Hermanson, J. C.; Choi, C. K.; Allen, J. S. Modeling Liquid–Vapor Phase Change Experiments: Cryogenic Hydrogen and Methane. Colloids and Surfaces A: Physicochemical and Engineering Aspects 2023, 675, 131932. https://doi.org/10.1016/j.colsurfa.2023.131932.
    (4) Bellur, K.; Médici, E. F.; Hussey, D. S.; Jacobson, D. L.; LaManna, J.; Leão, J. B.; Scherschligt, J.; Hermanson, J. C.; Choi, C. K.; Allen, J. S. Results from Neutron Imaging Phase Change Experiments with LH2 and LCH4. Cryogenics 2022, 125, 103517. https://doi.org/10.1016/j.cryogenics.2022.103517.


## Repository Contents

## Requirements

## Usage Instructions

## Citation and License