# Introduction
The scripts here calculates **Undulation Corrected** transverse electron density, mass density or neutron scattering length density (SLD) for lipid membranes or lipid membrane with different biomolecules like proteins or nucleic acids.  
## Why
1. Transverse electron or neutron SLD are one most important structural quantities for lipid bilayer. The scattering form factor obtained by fourier transform of these quantities or reflectivity profiles can be direclty compared with experiments to validate the forcefields predictions.
2. For Martini systems obtaining electron/mass density or neutron SLD profile calculation is non-trivial. First one needs to assign number of electrons to each bead, however, even after doing so using tools like ``gmx density`` is not feasible since a given bead maps to different number of electrons depending on the residue i.e lipid type or nucleotide or amino acid.
3. When simulating large membrane patches (>500 lipid per monolayer), the inherent undulations smear out the transverse density profiles. To make direct comparison with experiments those undulations need to be corrected.
4. The scripts here takes care of all the above for both **CG-Martini** and **All-atoms** lipid bilayer systems with and without other biomolecules.
   
