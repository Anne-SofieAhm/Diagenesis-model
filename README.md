# Diagenesis-model
Numerical model on fluid-rock interaction in carbonate rocks (Ahm et al. 2018)

Building on published modeling efforts (Banner & Hanson, 1990; Blattler & Higgins, 2015; Fantle & Higgins, 2014), we here present a MatLab based (R2016b) numerical model that simulates early marine diagenesis of carbonate sediments by stoichiometrically recrystallizing aragonite to low-Mg calcite (neomorphism) or aragonite/calcite to dolomite (dolomitization) assuming conservation of mass of carbon in the sediment. For a detailed description of the model setup we refer to Ahm et al. 2018 (GCA Vol. 236, pp 140-159, doi:10.1016/j.gca.2018.02.042).

To run the model you will need the following files:

(1) *Reactive_model_v6_base.m*

This is the main script used to run the model. In this script you define the initial conditions such as primary       mineralogy, diagentic mineralogy (m), reaction rate (R), flow rate (u), length scale of flow path (number of boxes), and composition of the diagenetic fluid. 

(2) *constants.m*

A script with some basic constants e.g. volume of each box and stoichiometry of the primary minerals. This script is also were you can change the composition of the primary mineral e.g. from aragonite to high-Mg calcite. 

(3) *dJ6dt.m*

The main ODE solver than calculates the change in mass and isotopic composition over time.

(4) *fluxes6.m*

This script calculates the input and output fluxes into each box at every time step. These fluxes are then pulled into the ODE solver (dJ6dt.m) to solve dM/dt and dM_delta_/dt.
