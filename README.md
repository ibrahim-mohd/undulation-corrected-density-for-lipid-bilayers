# Introduction
The scripts here calculates transverse electron density, mass density or neutron scattering length density (SLD) for lipid membranes or lipid membrane with different biomolecules like proteins or nucleic acids.

## How is mass or electron assigned to CG-Martini beads
The script `martini_bead_mass_electron_assigner.py` reads the Martini mapping file for a given molecule. A section of mapping file for `DOPC` is shown in the following Table
| **Index** | **Atoms** | **CG-Beads**               |
|-----------|-----------|----------------------------|
| 1         | N         | NC3                        |
| 2         | C12       | NC3 NC3 NC3 PO4            |
| 3         | C13       | NC3                        |
| 4         | C14       | NC3                        |
| 5         | C15       | NC3                        |
| 6         | H12A      | NC3 NC3 NC3 PO4            |
| 7         | H12B      | NC3 NC3 NC3 PO4            |
| …         | …         | …                          |

> **Table:** An example Martini mapping file. This portion is from the DOPC mapping file. Source: [Martini Mapping Files](https://cgmartini.nl/docs/downloads/force-field-parameters/martini2/lipidome.html)

---

Based on the mapping file (Table above), the $b_i$ for CG-bead **NC3** is assigned as:

$$
b_i (\mathrm{NC3}) = b_i (\mathrm{N}) + \frac{3 b_i (\mathrm{C})}{4} + 3b_i (\mathrm{C}) + \frac{3b_i (\mathrm{H})}{4} + \frac{3b_i (\mathrm{H})}{4} + \dots
$$

Where $b_i (X)$ is the scattering factor of atom $X$.
