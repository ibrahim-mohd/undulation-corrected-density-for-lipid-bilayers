# 1. Introduction
The scripts here calculate **Undulation Corrected** transverse **electron density**, **mass density**, or **neutron scattering length density (SLD)** for lipid membranes or lipid membranes with different biomolecules like proteins or nucleic acids.  

## 1.1 Why

&nbsp;&nbsp;&nbsp;&nbsp;1. Transverse electron or neutron SLD are among the most important structural quantities for lipid bilayers. The scattering form factor obtained by Fourier transform of these quantities or reflectivity profiles can be directly compared with experiments to validate the forcefield predictions.  

&nbsp;&nbsp;&nbsp;&nbsp;2. For Martini systems, obtaining electron/mass density or neutron SLD profile calculation is non-trivial. First, one needs to assign the number of electrons to each bead. However, even after doing so, using tools like `gmx density` is not feasible since a given bead maps to a different number of electrons depending on the residue (i.e., lipid type, nucleotide, or amino acid).  

&nbsp;&nbsp;&nbsp;&nbsp;3. When simulating large membrane patches (>500 lipids per monolayer), the inherent undulations smear out the transverse density profiles. To make direct comparisons with experiments, those undulations need to be corrected.  

&nbsp;&nbsp;&nbsp;&nbsp;4. The scripts here take care of all the above for both **CG‑Martini** and **All‑atom** lipid bilayer systems, with and without other biomolecules.  
## Dependencies
MDAnalysis (v-2.9.0), periodictable (v-2.0.2). Other version may work too
# 2. Martini Coarse-grained Membrane Simulations: undulation corrected density
There are two scripts that need to be run one after the other:  

1. **Assign electron/mass/neutron** to Martini beads using  `martini_bead_electron_mass_assigner.py`.

```bash
python martini_bead_electron_mass_assigner.py -m ./Mapping -c gro_file -o cg_bead_prop.json -dens electron
```

The `Mapping` directory is part of this repository obtained from the Martini website (backward tutorial zip file). ```-dens``` can be ``electron``, ``mass`` or ``neutron``
 

2. Obtain **undulation corrected/regular density** using `martini_membrane_uc_density.py`
```bash
python martini_membrane_uc_density.py -f mol.xtc -s npt.tpr -j cg_electron.json -uc 1 -o output_file
```
  The number of Fourier terms (-N), filter threshold (-q0), Group (-group) to consider for udulating reference surface is set to default values. You can find optimum values for your system by trying around or follow the paper for a procedure. 

# 3. All-atom Membrane Simulations: undulation corrected density
For all atoms system  to obtain **undulation corrected/regular density** only the `all_atom_membrane_uc_density.py` as:
```bash
python all_atom_membrane_uc_density.py -f mol.xtc -s npt.tpr -uc 1 -o output_file
```
  The number of Fourier terms (-N), filter threshold (-q0), Group (-group) to consider for udulating reference surface is set to default values. The default number of fourier term is N=4 which acutally means 4*2=8 terms. More than 4 terms just makes the code slower without any significant improvement in the calculations.

# 4. Determine the threshold wavevector ($q_0$) 
