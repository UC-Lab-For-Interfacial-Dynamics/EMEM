# Embedded Multiscale Evaporation Model (EMEM)

## Introduction
EMEM is a subgrid evaporation model for multiphase CFD simulations that provides the following:
- kinetic theory-based evaporation mass and heat calculation
- in situ calculation of the accommodation coefficient using transition state theory
- a lubrication theory-based transition thin film model to capture evaporation into the nanoscale region
- two-way coupling between the subgrid transition thin film model and the CFD simulation
- implementation within commercial-grade CFD software — Ansys Fluent^TM^.

## Theoretical Background
The development of this model is motivated by the multiscale nature of the menisci of wetting fluids, where significant evaporation flux is contributed by the transition thin film region and cannot be resolved within macroscale CFD simulations. EMEM employs a differential form of the augmented Young-Laplace equation to derive an ODE that describes the evolution of the interface shape. The ODE is solved numerically using a fixed-step fourth-order Runge-Kutta method at each iteration of the CFD simulation. The Navier-Stokes equations are reduced under the lubrication approximation and solved to obtain the hydrodynamic pressure. One-dimensional heat conduction from the solid wall to the interface is coupled with a latent heat loss condition based on kinetic evaporation. A two-way coupling is established with the CFD simulation:<br/>
1. Solid wall temperatures are extracted and used as boundary conditions at each CFD iteration
2. Source terms are returned to the CFD simulation to add the total evaporating mass as a source term at the average-interface temperature.

The code provided here is set up to simulate the extended meniscus of liquid Hydrogen in a 10 mm cylindrical container for validation with the experimental test case.

## References
See associated paper for details: <br/>
    [1] Yasin, A.; Pakanati, S.; Bellur, K. A multiscale CFD model of evaporating Hydrogen menisci: Incorporating subgrid thin-film dynamics and in situ accommodation coefficients. Fuels [Under Review]<br/><br/>
    [2] Yasin, A.; Bellur, K. Computational Modeling of Evaporation without Tuning Coefficients. Applied Thermal Engineering 2025, 276, 126807. https://doi.org/10.1016/j.applthermaleng.2025.126807.<br/><br/>
    [3] Bellur, K.; Médici, E. F.; Hermanson, J. C.; Choi, C. K.; Allen, J. S. Modeling Liquid–Vapor Phase Change Experiments: Cryogenic Hydrogen and Methane. Colloids and Surfaces A: Physicochemical and Engineering Aspects 2023, 675, 131932. https://doi.org/10.1016/j.colsurfa.2023.131932. <br/><br/>
    [4] Bellur, K.; Médici, E. F.; Hussey, D. S.; Jacobson, D. L.; LaManna, J.; Leão, J. B.; Scherschligt, J.; Hermanson, J. C.; Choi, C. K.; Allen, J. S. Results from Neutron Imaging Phase Change Experiments with LH2 and LCH4. Cryogenics 2022, 125, 103517. https://doi.org/10.1016/j.cryogenics.2022.103517. <br/><br/>
	[5] Ansys Fluent 12.0/12.1 Documentation (UDF Guide). [https://www.afs.enea.it/project/neptunius/docs/fluent/html/udf/main_pre.htm](https://www.afs.enea.it/project/neptunius/docs/fluent/html/udf/main_pre.htm).  <br/><br/>

## Repository Contents
- ThinFilm_MassSource.c
- ThinFilm_HeatBal.c

## Requirements
- Ansys Fluent^TM^ with User-Defined Functions (Code was developed and tested on v2023 R2).
- Built-in C compiler (Clang)

## Usage Instructions
1. Set up the Ansys Fluent mesh and simulation parameters with a *cutoff* at the extended meniscus region (see Ref [1]).
2. Configure the `ThinFilm_MassSource.c` and `ThinFilm_HeatBal.c` — Follow instruction in the file headers.
3. Compile the UDFs into the Ansys Fluent simulation file. (See Ansys Fluent UDF Guide).
4. Hook up the UDFs to the *cutoff* vapor mesh cell (see Ref[1]).
5. Execute the simulation. 

## Citation and License
EMEM is open-sourced under the [MIT License]([url](https://github.com/UC-Lab-For-Interfacial-Dynamics/EMEM/blob/main/LICENSE)).
Please cite our main paper (Ref. [1]) if you use EMEM in your work.

## Development
Please report bugs to Ayaaz Yasin ([yasinaz@mail.uc.edu](mailto:yasinaz@mail.uc.edu)).<br/>
Forthcoming updates: code for the non-tuning bulk meniscus evaporation model.


