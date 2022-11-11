The existing LRR Reynolds stress model in OpenFOAM uses the simplified
rapid pressure strain model and linear wall damping function which 
previously produced weaker turbulence anisotropy and secondary currents.
Therefore, the full LRR pressure strain model, the nonlinear wall & 
free surface damping functions, and the free surface boundary condition 
for the rate of turbulent kinetic energy dissipation were incorporated
in OpenFOAM based on literature.

The modifications were tested in -Dev version. They should be working 
in OpenFOAM-10 also. 

For more details follow the manuscript titled 
"Reynolds Stress Modeling of Supercritical Narrow Channel Flows using OpenFOAM: Secondary Currents and Turbulent Flow Characteristics"
accepted in Physics of Fluids journal. https://doi.org/10.1063/5.0124076
