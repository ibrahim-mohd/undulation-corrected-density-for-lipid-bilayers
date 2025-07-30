#author: Mohd Ibrahim, Technical University of Munich
# email: ibrahim.mohd@tum.de

import argparse
import json
import logging
from pathlib import Path
import numpy as np
import MDAnalysis as mda
import periodictable as pt
import warnings

# Suppress specific warnings from MDAnalysis
warnings.filterwarnings("ignore")

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# Default element data
ELEMENTS =  {
    # Most atoms in protein/lipid and nucleic acids
    # if required you can simply add more atoms here
    "H":  {"name": "Hydrogen",   "electron":pt.H.number,  "mass": pt.H.mass, "neutron":pt.H.neutron.b_c},
    "D":  {"name": "Deuterium",  "electron":pt.D.number,  "mass": pt.D.mass, "neutron":pt.D.neutron.b_c},
    "C":  {"name": "Carbon",     "electron":pt.C.number,  "mass": pt.C.mass, "neutron":pt.C.neutron.b_c},
    "N":  {"name": "Nitrogen",   "electron":pt.N.number,  "mass": pt.N.mass, "neutron":pt.N.neutron.b_c},
    "O":  {"name": "Oxygen",     "electron":pt.O.number,  "mass": pt.O.mass, "neutron":pt.O.neutron.b_c},
    "P":  {"name": "Phosphorus", "electron":pt.P.number,  "mass": pt.P.mass, "neutron":pt.P.neutron.b_c},
    "S":  {"name": "Sulfur",     "electron":pt.S.number,  "mass": pt.S.mass, "neutron":pt.S.neutron.b_c},
    # Common ions 
    "Na": {"name": "Sodium",     "electron":pt.Na.number,  "mass": pt.Na.mass, "neutron":pt.Na.neutron.b_c},
    "K":  {"name": "Potassium",  "electron":pt.K.number,  "mass": pt.K.mass, "neutron":pt.K.neutron.b_c},
    "Ca": {"name": "Calcium",    "electron":pt.Ca.number,  "mass": pt.Ca.mass, "neutron":pt.Ca.neutron.b_c},
    "Mg": {"name": "Magnesium",  "electron":pt.Mg.number,  "mass": pt.Mg.mass, "neutron":pt.Mg.neutron.b_c},
    "Cl": {"name": "Chlorine",   "electron":pt.Cl.number,  "mass": pt.Cl.mass, "neutron":pt.Cl.neutron.b_c}
}



def load_mapping_file(mapping_dir: Path, resname: str) -> list[str]:
    """Try to load mapping file for a given residue name."""
    candidates = [
        mapping_dir / f"{resname.lower()}.charmm36.map",
        mapping_dir / f"{resname.lower()}.amber.map"
    ]
    for candidate in candidates:
        if candidate.exists():
            logging.info(f"Using mapping file: {candidate}")
            with open(candidate, "r") as f:
                return f.readlines()
    raise FileNotFoundError(
        f"No mapping file found for {resname}. "
        f"Tried: {candidates[0].name}, {candidates[1].name}"
    )


def build_default_beads(dens_type: str) -> dict:
    """Define default beads for water and ions."""
    H, O = ELEMENTS["H"][dens_type], ELEMENTS["O"][dens_type]
    return dict(
        W=4 * (O + 2 * H),
        WF=4 * (O + 2 * H),
        NA=ELEMENTS["Na"][dens_type] + 6 * (O + 2 * H),
        CL=ELEMENTS["Cl"][dens_type] + 6 * (O + 2 * H),
        MG=ELEMENTS["Mg"][dens_type],
        K=ELEMENTS["K"][dens_type],
        CA=ELEMENTS["Ca"][dens_type],
        WP=0,
        WM=0,
    )


def parse_mapping_file(mapping_lines: list[str], dens_type: str) -> dict:
    """Parse mapping file lines and return CG bead assignments."""
    start_reading = False
    mapping_list = []
    cg_beads = set()
    molecule_name = None

    for idx, line in enumerate(mapping_lines):
        if "[ molecule ]" in line:
            molecule_name = mapping_lines[idx + 1].split(";")[0].strip().split()[0]

        if "[ atoms ]" in line:
            start_reading = True
            continue

        if start_reading:
            if "[" in line:
                break
            stripped = line.split(";")[0].strip().split()
            try:
                int(stripped[0])  # ensure first element is atom index
                if len(stripped) > 2:
                    cg_beads.update(stripped[2:])
                mapping_list.append(stripped)
            except (ValueError, IndexError):
                pass

    bead_dict = {bead: 0 for bead in cg_beads}
 
    for line in mapping_list:
        
        try:
            
            n_electron = ELEMENTS[line[1][0]][dens_type]
            # Note that int he mapping file sometimes in a given line there is no CG bead
            # in those cases the atom belongs to the beads in the line before
            if len (line) >2:
                cg_bead_in_line = line [2:]
            
            electron_per_cg_bead = n_electron/len(cg_bead_in_line)
            
            for key in cg_bead_in_line: bead_dict [key] += electron_per_cg_bead
            
        
        except (KeyError): pass

    return molecule_name, bead_dict, mapping_list


def main():
    parser = argparse.ArgumentParser(
        description="Assign electron/mass/neutron scattering factors to Martini beads by reading mapping file"
    )
    parser.add_argument("-m", dest="mapping_file_path", type=str, default="./Mapping",
                        help="Martini mapping file directory")
    parser.add_argument("-c", dest="gro_file", type=str, default="npt.gro",
                        help="input structure file, gro or pdb")
    parser.add_argument("-o", dest="out_file", type=str, default="cg_electron.dat",
                        help="output file")
    parser.add_argument("-dens", dest="dens_type", type=str, default="electron",
                        choices=["electron", "mass", "neutron"], help="density type can be electron or mass or neutron coherent scattering lengths")
    args = parser.parse_args()

    mapping_path = Path(args.mapping_file_path)
    dens_type = args.dens_type

    u = mda.Universe(args.gro_file)

    # we select all non-solvnet resnames
    # the solvnet electron are assigned by default
    resnames = [
        r for r in set(u.atoms.residues.resnames)
        if r not in {"W", "WF", "PW", "ION", "CL", "NA"}
    ]

    cg_bead_electron = {res: {} for res in resnames}
    default_beads = build_default_beads(dens_type)

    for resname in resnames:
        try:
            mapping_lines = load_mapping_file(mapping_path, resname)
            molecule_name, bead_dict, mapping_list = parse_mapping_file(mapping_lines, dens_type)
            cg_bead_electron[molecule_name] = bead_dict

            total_val = int(np.round(sum(bead_dict.values())))
            logging.info(f"{molecule_name}: {len(mapping_list)} atoms, total {dens_type} = {total_val}")

        except FileNotFoundError as e:
            logging.warning(e)

    # Add default solvent/ion beads
    cg_bead_electron["SOL"] = default_beads

    # Write JSON output
    with open(Path(args.out_file).with_suffix(".json"), "w") as f:
        json.dump(cg_bead_electron, f, indent=4)
 
    logging.info(f"Processed {len(resnames)} non-solvent residues.")

    # make sure all the CG residues has been assigned electrons/mass
    
    all_resnames = set(u.atoms.residues.resnames)
    residues_assigned_electron = list(cg_bead_electron.keys()) + list(cg_bead_electron ["SOL"].keys())
    not_account_for = 0
    for res in all_resnames:
        if res not in residues_assigned_electron and res not in ["ION", "PW"]:
            print (f"Residue {res} not assigned {dens_type}")
            not_account_for += 1
    if not_account_for: print (f"{not_account_for} residues not assigned {dens_type}")
    else: print (f"All residues has been assigned {dens_type}")

if __name__ == "__main__":
    main()

