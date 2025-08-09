# 1. Introduction
The scripts here calculate **Undulation Corrected** transverse **electron density**, **mass density**, or **neutron scattering length density (SLD)** for lipid membranes or lipid membranes with different biomolecules like proteins or nucleic acids. The undulation correction procedure is performed by using the method from Braun et al.

*Braun, A.R., Brandt, E.G., Edholm, O., Nagle, J.F. and Sachs, J.N., 2011. Determination of electron density profiles and area from simulations of undulating membranes. Biophysical journal, 100(9), pp.2112-2120*.

## 1.1 Why

&nbsp;&nbsp;&nbsp;&nbsp;1. Transverse electron or neutron SLD are among the most important structural quantities for lipid bilayers. The scattering form factor obtained by Fourier transform of these quantities or reflectivity profiles can be directly compared with experiments to validate the forcefield predictions.  

&nbsp;&nbsp;&nbsp;&nbsp;2. For Martini systems, obtaining electron/mass density or neutron SLD profile calculation is non-trivial. First, one needs to assign the number of electrons to each bead. However, even after doing so, using tools like `gmx density` is not feasible since a given bead maps to a different number of electrons depending on the residue (i.e., lipid type, nucleotide, or amino acid).  

&nbsp;&nbsp;&nbsp;&nbsp;3. When simulating large membrane patches (>500 lipids per monolayer), the inherent undulations smear out the transverse density profiles. To make direct comparisons with experiments, those undulations need to be corrected.  

&nbsp;&nbsp;&nbsp;&nbsp;4. The scripts here take care of all the above for both **CG‑Martini** and **All‑atom** lipid bilayer systems, with and without other biomolecules.  
## 1.2 Dependencies
MDAnalysis (v-2.9.0), periodictable (v-2.0.2). Other version  should work as well since I am not using anything version specific
# 2. Martini Coarse-grained Membrane Simulations: undulation corrected density
There are two scripts that need to be run one after the other:  

1. **Assign electron/mass/neutron** to Martini beads using  `martini_bead_electron_mass_assigner.py`.

```bash
python martini_bead_electron_mass_assigner.py -m ./Mapping -c gro_file -o cg_bead_prop.json -dens electron
```

The `Mapping` directory is part of this repository obtained from the Martini website (backward tutorial zip file). ```-dens``` can be ``electron``, ``mass`` or ``neutron``
 

2. Obtain **undulation corrected/regular density** using `martini_membrane_uc_density.py`
```bash
python martini_membrane_uc_density.py -f mol.xtc -s npt.tpr -j cg_electron.json -uc 1 -o output_file -q0 0.04 -N 4 -group "C3A"
```
  The number of Fourier terms (-N), wave vector threshold (-q0), Group (-group) to consider for udulating reference surface is set to default values.

# 3. All-atom Membrane Simulations: undulation corrected density
For all atoms system  to obtain **undulation corrected/regular density** only the `all_atom_membrane_uc_density.py` needs to be run as:
```bash
python all_atom_membrane_uc_density.py -f mol.xtc -s npt.tpr -uc 1 -o output_file -q0 0.04 -N 4 -group "C114"
```
  The number of Fourier terms (-N), filter threshold (-q0), Group (-group) to consider for udulating reference surface is set to default values. The default number of fourier term is N=4 which acutally means 4*2=8 terms. More than 4 terms just makes the code slower without any significant improvement in the calculations.


<img src="/Figures/dppc_dopc_edens.png" width="1000"> 

Figure 1: <em> (A) Electron density of pure DPPC bilayer from CG-Martini simulations. There is no undulation smearing for a small bilayer with $N_{lip}=100$ lipids/monolayer. For $N_{lip}=1000$, undulations smear out the electron density. After undulation correction we recover the electron density that is closer to $N{lip}=100$ (B) These profiles are obtained from all-atom projected of Martini DOPC/tRNA system therefore usign the all-atom undulation correction script.</em> 


# 4. Determination of the threshold wavevector ($q_0$) 
For the undulation correction choosing proper $q_0$ is critcal. If its too high, we are correcting for thermal motion and if two low we omit the acutal inherant undulations. To obtain the optimum we follow the procedure by Braun et al., where we obtain the **spectral intensity** by scanning a relatively large $q-$ space i.e $N>=50$ or so. For membranes, the spectral intesity of actual unduation scales as $q^{-4}$, so we choose the $q_0$ where the behaviour deviates from this scaling. In most case of DOPC and DPPC for both all-atom and martini, the default value i found i.e $q_0 \sim 0.04~Å^{-1}$ works well. However, to decide on a given value, always **simulate a small bilayer (~64 lipids/monolayer)** for few tens of nanoseconds. The best $q_0$ should produce denisty profiles that is closest to the small system.  To explore the spectral intensity behaviour for different sytems, you can use the provided notebook `obtain_threshold_wave_vector.ipynb`.

<img src="/Figures/dppc_1000_spectrum.png" width="1000"> 

Figure 2: <em> Undulation spectra for pure DPPC bilayer from CG-Martini simulations with 1000 lipids/monolayer. At low frequency, the undulations scale as $q^{-4}$, these corresponds to actual undulations regimes. At high-$q$, the intensity starts increaing again, this regime has to be filtered out. (A) without filter (B) after applying filter from Braun et al. Here we used two different atom groups PO4 and C3A to calculate the undulation reference surface</em>

# Remarks
1. The density calculation output has always 4 columns: (i) z-distance (ii) Total density profile (iii) All not solvent and ions i.e, (lipids + else) and (iv) solvent + ions.
2. You do not need to specify the solvent: For all-atom systems the code assumes your solvent has one of the following residue names, taken from MDAnalysis documentation:
         ``` ["H2O", "HOH", "OH2", "HHO", "OHH", "TIP", "T3P", "T4P", "T5P", "SOL", "WAT", "TIP2", "TIP3", "TIP4"]```
 
   
And ions names is one of the following: ```["NA", "CL", "K","CA", "SOD","POT","CAL", "CLA","Cl-", "Na+","K+", "Li+", "MG", "NIO", "CXY", "CIO", "LIO", "KIO", "mMg", "nMg"]```

4. If you get error related to sol/ion resname not found, just add the names of your ion or solvnet to the above list in the script.

# 5. References
If you find this thing useful please cite our work: 

Mohd Ibrahim, Jürgen Kofinger, Martin Zacharias, Emanuel Schneck, and Nadine Schwierz "RNA at lipid/water interfaces: Probing interactions and structures" (*In preparation*)
