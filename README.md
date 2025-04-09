# N-Body Simulations for Planetary Systems and Arenstorf Translunar Insertions

This project describes the implementation and analysis of the Forward Euler (FE), Heun, Kick-Drift-Kick (KDK), and RK4 numerical methods to simulate various N-body problems with the goal of being able to solve for Arenstorf orbits to enable a Translunar Insertion (TLI). After analyzing the stability regions of these numerical methods, examining the effect that time step size and truncation error have on various simulations, and validating their accuracy by comparing their results to real planetary orbit data from JPL, the stable TLI used by the Apollo missions was solved for and replicated accurately.

Additionally, this project includes code for:

* **RK4 Moon Orbit Estimation:** Simulates the Earth-Moon system using the RK4 method and compares the estimated Moon orbit with real data obtained from JPL's ephemeris. This includes handling proper rotation for accurate initial conditions and visualizing the 3D and planar orbits.
* **Stability Region Analysis:** Generates plots illustrating the stability regions of the Forward Euler, Heun, RK4, and Leapfrog numerical methods in the complex plane.
* **Two-Body Simulation:** Implements an RK4 simulation for a general two-body system, allowing for the analysis of the effect of varying time step sizes on the orbital dynamics.

## Code

The project code is implemented in MATLAB. It includes several scripts for different simulation scenarios and analyses:

* **Main TLI Simulation Script:** Contains the implementation of the Forward Euler, Heun, KDK, and RK4 methods for solving the Arenstorf TLI problem.
* **Four-Body Planetary System Simulation Script:** Simulates the orbits of the Sun, Mercury, Earth, and Mars using the four numerical methods.
* **`RK4_Moon_Orbit_Estimation.m`:** Simulates the Earth-Moon system using the RK4 method.
* **`Stability_region.m`:** Generates stability plots for the numerical methods.
* **`Two_body_simulation.m`:** Simulates a general two-body system with varying time step sizes.

## Functions

The MATLAB code includes the following key functions:

* `du_dt(x, x_t, y, y_, mu, mu_t)`: Calculates the acceleration in the x-direction for the circular restricted three-body problem.
* `dv_dt(x, x_, y, y_, mu, mu_t)`: Calculates the acceleration in the y-direction for the circular restricted three-body problem.
* `dx_dt(u)`: Returns the velocity in the x-direction.
* `dy_dt(v)`: Returns the velocity in the y-direction.
* `du(u)`: Defines the system of first-order ordinary differential equations for the four-body problem (Sun, Mercury, Earth, Mars).
* `du_KDK(X)`: Defines the acceleration for each of the four bodies in the KDK method.
* **`dv_dt(m_j, r_i, r_j)` (in `RK4_Moon_Orbit_Estimation.m`):** Calculates the gravitational acceleration between two bodies.
* **`dx_dt(m_j, r_i, r_j, v_i)` (in `RK4_Moon_Orbit_Estimation.m`):** Returns the velocity of a body.
* **`dV_dT(m_j, r_i, r_j)` (in `Two_body_simulation.m`):** Calculates the gravitational acceleration between two bodies (2D).
* **`dX_dT(m_j, r_i, r_j, v_i)` (in `Two_body_simulation.m`):** Returns the velocity of a body (2D).

## Simulation Examples

The code simulates the following scenarios:

1.  **Four-Body Planetary System (Sun, Mercury, Earth, Mars):** This simulation uses real initial position and velocity data obtained from JPL's ephemeris to model the orbits of these inner solar system bodies. The results from the four numerical methods are compared visually and quantitatively (by calculating the percentage difference in the orbital radius of Mercury compared to the "true" orbit obtained from JPL data).

2.  **Arenstorf Translunar Insertion (TLI):** This simulation focuses on replicating the stable TLI orbit used by the Apollo missions within the circular restricted three-body problem (Earth-Moon system). The code solves for a specific initial condition that leads to a successful lunar insertion trajectory. The results are visualized, showing the probe's trajectory relative to the Earth and Moon.

3.  **Earth-Moon System (RK4):** This simulation uses the RK4 method to model the orbit of the Moon around the Earth, utilizing real initial conditions and allowing for comparison with actual lunar trajectory data.

4.  **General Two-Body System (RK4):** This simplified simulation models the interaction of two bodies in a plane, enabling the study of how the time step size affects the accuracy and stability of the numerical solution.

## Plots

The MATLAB scripts generate several figures:

* **TLI.png:** Visualizes the Arenstorf TLI orbit obtained using the RK4 method, showing the trajectory relative to the Earth and Moon.
* **Figure 1:** Compares the simulated orbits of Mercury, Earth, and Mars around the Sun using the Forward Euler, Heun, and RK4 methods.
* **Figure 46:** Plots the percentage difference in the orbital radius of Mercury over time for the Heun, KDK, and RK4 methods compared to the JPL data.
* **Figure 7:** Shows the Arenstorf TLI orbit obtained using the Forward Euler method.
* **Stability.jpg (Figure 1 in `Stability_region.m`):** Illustrates the stability regions of the Forward Euler, Heun, RK4, and Leapfrog numerical methods in the complex plane.
* **Figure 1 (in `RK4_Moon_Orbit_Estimation.m`):** Shows the 3D simulated orbit of the Moon around the Earth compared to the true orbit data.
* **Figure 2 (in `RK4_Moon_Orbit_Estimation.m`):** Shows the planar (2D) simulated orbit of the Moon around the Earth compared to the true orbit data.
* **N=2\_Simulation.png (Figure 1568 in `Two_body_simulation.m`):** Displays the orbits of two bodies for different time step sizes.

## Getting Started

To run the simulations, you will need:

1.  **MATLAB:** Ensure you have MATLAB installed on your system. The code utilizes the `planetEphemeris` function, which requires the Aerospace Toolbox.
2.  **Aerospace Toolbox (Optional but Recommended):** While the core numerical methods are implemented from scratch, the `planetEphemeris` function from this toolbox is used to obtain real planetary and Earth-Moon orbit data for validation. If you do not have this toolbox, you would need to provide your own initial conditions for the astronomical simulations.
3.  **Download the Code:** Clone or download the MATLAB scripts (`.m` files) from this repository.
4.  **Run the Scripts:** Open each MATLAB script in the MATLAB environment and run them. Each script will execute its specific simulation and generate the corresponding plots.

**Note:** The initial conditions and simulation parameters (time steps, total time, number of steps) are defined at the beginning of each script and can be adjusted to explore different scenarios or improve accuracy. The plotting commands will automatically generate the figures. For the `RK4_Moon_Orbit_Estimation.m` script, the rotation calculations aim to align the initial conditions with the orbital plane.
