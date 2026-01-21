# Overview
Though many one-dimensional models exist, the one-dimensional turbulence (ODT; Wunch and Kerstein, 2005) model has been demonstrated to simulate the dry Rayleigh Bernard convection (RBC) for very high Rayleigh numbers, compared to lab experiments. The model parameters were tuned to match the laboratory observations of the Nusselt number, Sherwood number, and scaling of these non-dimensional quantities in high Rayleigh number experiments.

The ODT model is a stochastic turbulence model that randomly selects an eddy size from an eddy-size distribution determined by the RBC boundary conditions. The effect of an eddy is implemented through triplet mapping and numerical diffusion. Triplet mapping rearranges grid cells within an eddy of length l, such that the scalar field has a sinusoidal variation. This increases the gradients between adjacent cells and increases the scalar diffusion. For RBC convection under steady-state conditions, the main variability is along the vertical direction. Therefore, the one-dimensional model is set up along the vertical direction.

This model simulates dry Rayleigh-Bénard convection with only temperature as its scalar field. To simulate temperature and water vapor, use moist-ODT whereas to simulate temperature, water vapor, and cloud droplets, use cloudy-ODT.

# Contributions
* **Scott Wunsch and Alan Kerstein** developed the dry ODT model. 
* **Mani Rajagopal** modified the code to output eddy information, scalar probability distribution, and evolution of domain statistics. Also included is code for input and output to the respective folders, as well as a random initial perturbation seed (as command line argument) to run ensemble of experiments.

# References
Wunsch, S., & Kerstein, A. R. (2005). A stochastic model for high-Rayleigh-number convection. Journal of Fluid Mechanics, 528, 173-205.
