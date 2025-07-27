# Written by Mohd Ibrahim
# Technical University of Munich
# Email: ibrahim.mohd@tum.de
import numpy as np
import json
import argparse

parser = argparse.ArgumentParser(description="Assign electron/mass to Martini beads by reading mapping file")
parser.add_argument('-m', dest='mapping_file', type=str, default='./Mapping/dopc.charmm36.map',help='Martini mapping file')
parser.add_argument('-o', dest='out_file', type=str, default='cg_electron.dat',help='output file')
parser.add_argument('-dens', dest='dens_type', type=str, default='electron',help='density type electron or mass')
  
args      = parser.parse_args()
in_file   = args.mapping_file
out_file  = args.out_file
dens_type = args.dens_type
#### density either mass or electron
assert dens_type in ["electron", "mass"]

element_list = {
    # Most atoms in protein/lipid and nucleic acids
    # if required you can simply add more atoms here
    "H":  {"name": "Hydrogen",    "electron": 1,  "mass": 1.008},
    "C":  {"name": "Carbon",      "electron": 6,  "mass": 12.011},
    "N":  {"name": "Nitrogen",    "electron": 7,  "mass": 14.007},
    "O":  {"name": "Oxygen",      "electron": 8,  "mass": 15.999},
    "P":  {"name": "Phosphorus",  "electron": 15, "mass": 30.974},
    "S":  {"name": "Sulfur",      "electron": 16, "mass": 32.06},
    # Common ions 
    "Na": {"name": "Sodium",      "electron": 11, "mass": 22.990},
    "K":  {"name": "Potassium",   "electron": 19, "mass": 39.098},
    "Ca": {"name": "Calcium",     "electron": 20, "mass": 40.078},
    "Mg": {"name": "Magnesium",   "electron": 12, "mass": 24.305},
    "Cl": {"name": "Chlorine",    "electron": 17, "mass": 35.45}
}

# dictionary of default beads like water and ions
# 4 water molecules map to 1 bead
# For ions, ions are modelled with their hydration shell. For now the hydration number for Na and Cl is considered. 6 for both
#

hydrogen   =  element_list["H"][dens_type]
oxygen     =  element_list["O"][dens_type]
sodium     =  element_list["Na"][dens_type]
chlorine   =  element_list["Cl"][dens_type]
magnesium  =  element_list["Mg"][dens_type]
potassium  =  element_list["K"][dens_type]
calcium    =  element_list["Ca"][dens_type]


default_beads = dict (W=4*(oxygen+2*hydrogen), 
                           WF=4*(oxygen+2*hydrogen),
                           NA = sodium + 6*(oxygen+2*oxygen),
                           CL = chlorine + 6*(oxygen+2*oxygen),
                           MG = magnesium,
                           K  = potassium,
                           CA = calcium)

################################################
mapping_file = open (in_file, 'r+').readlines ()

start_reading   = False
CG_beads_names  = set ()
mapping_list    = [] 

for index, line in enumerate (mapping_file):
    
    # get the CG molecule name for sanity check 
    if "[ molecule ]" in line: molecule_name =  mapping_file [index +1].strip()
    
    # the all-atom to CG mappign starts after "[ atoms ]" line
    if "[ atoms ]" in line:
        start_reading = True
        continue
        
    if start_reading:
        
        if "[" in line: break # Here the next section with other information starts so we stop here

        stripped = line.split (";") [0].strip ().split ()
        
        try:
            atom_index = int (stripped[0]) # The first element is an integer and is the index of atom
             
            # There are lines with hydrogen for which the CG bead is the explicityly mentioned
            # In theose cases the hydrogen mapps to CG beads in the preceding line
            
            if len (stripped) > 2:
                for name in stripped [2:]: CG_beads_names.add (name)
                
            mapping_list.append (stripped)
            
        except (ValueError, IndexError):   pass
            
# Dictionary for CG-beads electron
Bead_electron = {molecule_name:{keys:0 for keys in CG_beads_names}}

for line in mapping_list:
    try:
        
        n_electron = element_list [line[1][0]] [dens_type]
        
        if len (line) >2:
            cg_bead_in_line = line [2:]
        
        electron_per_cg_bead = n_electron/len(cg_bead_in_line)
        
        for key in cg_bead_in_line: Bead_electron [molecule_name][key] += electron_per_cg_bead
        
        
    except (KeyError):
        
        pass

total_number_of_electron_assigned = int (np.round (sum (Bead_electron[molecule_name].values ())))

## Write a json output file
Bead_electron ["SOL"] = default_beads
with open(out_file.split (".")[0]+".json", "w") as f:
    json.dump(Bead_electron, f, indent=4)

    
# write .dat output as input for groamcs "gmx density function"
number_of_cg_beads =  len(Bead_electron [molecule_name]) + len (default_beads)

with open(out_file, "w") as f:
    f.write (f"{number_of_cg_beads}\n")
    
    for key, value in Bead_electron[molecule_name].items(): f.write (f"{key} = {value}\n")

    for key, value in default_beads.items(): f.write (f"{key}={value}\n")
   
     
print (f"Total number of atoms in {molecule_name}: {len(mapping_list)}")
print (f"Total {dens_type} assigned to CG {molecule_name}: {total_number_of_electron_assigned}")
print (f"Please check if the total {dens_type} of/in a {molecule_name} molecule = {total_number_of_electron_assigned}")


