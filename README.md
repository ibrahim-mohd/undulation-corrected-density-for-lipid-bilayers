# 1. Introduction
The scripts here calculate **Undulation Corrected** transverse **electron density**, **mass density**, or **neutron scattering length density (SLD)** for lipid membranes or lipid membranes with different biomolecules like proteins or nucleic acids.  

## 1.1 Why

&nbsp;&nbsp;&nbsp;&nbsp;1. Transverse electron or neutron SLD are among the most important structural quantities for lipid bilayers. The scattering form factor obtained by Fourier transform of these quantities or reflectivity profiles can be directly compared with experiments to validate the forcefield predictions.  

&nbsp;&nbsp;&nbsp;&nbsp;2. For Martini systems, obtaining electron/mass density or neutron SLD profile calculation is non-trivial. First, one needs to assign the number of electrons to each bead. However, even after doing so, using tools like `gmx density` is not feasible since a given bead maps to a different number of electrons depending on the residue (i.e., lipid type, nucleotide, or amino acid).  

&nbsp;&nbsp;&nbsp;&nbsp;3. When simulating large membrane patches (>500 lipids per monolayer), the inherent undulations smear out the transverse density profiles. To make direct comparisons with experiments, those undulations need to be corrected.  

&nbsp;&nbsp;&nbsp;&nbsp;4. The scripts here take care of all the above for both **CG‑Martini** and **All‑atom** lipid bilayer systems, with and without other biomolecules.  

# 2. Martini Coarse-grained Membrane
There are two scripts that need to be run one after the other:  

1. The script `martini_bead_electron_mass_assigner.py` reads a `.gro` or `.pdb` file of the full system and, using the `Mapping` directory, assigns electron/mass/neutron scattering factors for each bead and creates a `.json` file.

```bash
       python martini_bead_electron_mass_assigner.py -m ./Mapping -c gro_file -o cg_bead_prop.json -dens electron
 
