#author: Mohd Ibrahim, Technical University of Munich
# email: ibrahim.mohd@tum.de
import argparse
import json
import pickle
import warnings

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element

import numpy as np
import periodictable as pt
from tqdm import tqdm

# Configure warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser (description="This code calculates transverse mass or electron or neutron scattering lenght density profiles for membrane or membrane/protein/nucleic acids system from all-atom MD simulations. It also perform undulation correction for comparison with expeirments especially ")

################  ###########################################

parser.add_argument('-f', dest='xtc_file', type=str, default='all.xtc',help='xtc file')
parser.add_argument('-s', dest='tpr_file', type=str, default='em.tpr',help='input tpr file')
parser.add_argument('-uc', dest='undulation_correct', type=int, default=0, help='if set to 1, undulation correction is performed')
parser.add_argument("-dens", dest="dens_type", type=str, default="electron",
                        choices=["electron", "mass", "neutron"], help="density type can be electron or mass or neutron coherent scattering lengths")
    
parser.add_argument('-dz', dest='dz', type=float, default= 0.25, help='slice width in Angstrom')

parser.add_argument('-q0', dest='q0', type=float, default= 0.04, help='filter boundary [A^-1]')

parser.add_argument('-N', dest='N', type=int, default= 4, help='number of fourier terms')

parser.add_argument('-skip', dest='skip', type=int, default= 1, help='skip frame by this interval')

parser.add_argument('-b', dest='begin_time', type=int, default= 0, help='Time of first frame to read from trajectory (unit ps)')
parser.add_argument('-e', dest='end_time', type=int, default= None, help='Time of last frame to read from trajectory (unit ps)')

parser.add_argument('-o', dest='output_filename', type=str, default=None,help='output filename')

parser.add_argument('-group', dest='group', type=str, default='C114',help='atoms to  consider for fourier coefficient calcualtion')
parser.add_argument('-z_lim', dest='z_lim', type=float, default= None,help='z limit: between -z_lim to +z_lim')
 
        
def get_fourier_coeff(box_dimensions, R_upper, R_lower, N=4, M=4, q0=1.15/10, filter=True):
    Lx, Ly = box_dimensions [0], box_dimensions [1]

    # Generate Q vectors (qx, qy) and their magnitudes
    n_vals = np.arange(-N, N+1)
    m_vals = np.arange(-M, M+1)
    n, m = np.meshgrid(n_vals, m_vals)
    n, m = n.flatten(), m.flatten()
    
    Q = np.column_stack((2 * n * np.pi / Lx, 2 * m * np.pi / Ly))
    Q = Q[np.any(Q != 0, axis=1)]  # Remove (0,0)
    
    Q_mag = np.sqrt(Q[:, 0]**2 + Q[:, 1]**2)
    
    # Precompute G
    if filter:
        G = np.sqrt(1 / (1 + (Q_mag / q0)**4))
    else: G = np.ones (len(Q_mag))
    # Combine R_upper and R_lower
    R_combined = np.vstack((R_upper, R_lower))
    
    # Extract x, y, and z components
    r_xy = R_combined[:, :2]  # x and y components
    z_combined = R_combined[:, 2]  # z component
    
    # Compute Q dot r for each Q and each r in R_combined
    Q_dot_r = np.dot(Q, r_xy.T)
    
    # Calculate the real and imaginary parts of the Fourier coefficients
    cos_QR = np.cos(Q_dot_r)
    sin_QR = np.sin(Q_dot_r)
    
    U_real = np.dot(cos_QR, z_combined)
    U_imag = -np.dot(sin_QR, z_combined)
    
    # Normalize
    norm_factor = len(R_combined)
    U_real /= norm_factor
    U_imag /= norm_factor
    
    # Combine the real and imaginary parts into complex numbers
    U = U_real + 1j * U_imag
    
    # Apply the filter
    U_filtered = G * U
    
    return Q, U_filtered


def inverse_transform(r, Q, U_filtered):
    
    Q = np.array (Q)
    # Calculate the dot product of r with each q in Q (vectorized)
    dots = np.dot(Q, r)

    # Calculate the exponential terms (vectorized)
    exp_terms = np.exp(1j * dots)

    # Compute u_tilde (vectorized)
    u_tilde = np.sum(U_filtered * exp_terms)

    # Compute delta_u_tildex and delta_u_tildey (vectorized)
    qx = Q[:, 0]
    qy = Q[:, 1]
    
    delta_u_tildex = np.sum(1j * qx * U_filtered * exp_terms)
    delta_u_tildey = np.sum(1j * qy * U_filtered * exp_terms)

    delta_u = delta_u_tildex + delta_u_tildey

    cos_theta = 1 / np.sqrt(delta_u.imag**2 + delta_u.real**2 + 1)

    return u_tilde, cos_theta

def get_different_monolayers (u, group):
    # To divide resides into bilayer leaflets, we first find the residues to which the group picked
    # for undulating reference surface belongs. Then we find the atoms with hightes and lowest 
    # z-coordinate. We use those atoms to seperate the leaftlet
    # to be sure we try both the highest z and lowest z and later on pick the atom which gives
    # most symmetric difference between the two. In most cases bothh should give same results
    res_ids          = " ".join (str (x) for x in u.select_atoms (f"name {group}").residues.resids)
    sel_top          = u.select_atoms(f"resid {res_ids} and not name H* 1H* 2H*")
    
    atom_name_min    =  sel_top.names [np.argmin (sel_top.positions [:,2])]
    atom_name_max    =  sel_top.names [np.argmax (sel_top.positions [:,2])]
    
    difference = {atom_name_min:0, atom_name_max:0}
    
    for name in [atom_name_min, atom_name_max]:
        
        top_bottom_layer = u.select_atoms (f"name {name}")
        z_com            =  top_bottom_layer.center_of_mass () [2]
        
        upper_layer = np.unique(u.select_atoms (f"name {name} and prop z>{z_com}").residues.resids)
        lower_layer = np.unique(u.select_atoms (f"name {name} and prop z<{z_com}").residues.resids)
        
        difference [name] = np.abs(len (upper_layer)-len(lower_layer))
        
    ################### select the criteria (highest or lowest) based on how symemtric the divide is
    # we assume a possible more symmetirc bilayer; but should not matter for most of the systems
    name_leaflet_boundary = min(difference, key=difference.get) 
    top_bottom_layer = u.select_atoms (f"name {name_leaflet_boundary}")
    z_com            =  top_bottom_layer.center_of_mass () [2]
    
    upper_layer = np.unique(u.select_atoms (f"name {name_leaflet_boundary} and prop z>{z_com}").residues.resids)
    lower_layer = np.unique(u.select_atoms (f"name {name_leaflet_boundary} and prop z<{z_com}").residues.resids)
    
    upper_layer_resid = " ".join (str(x) for x in upper_layer)
    lower_layer_resid = " ".join (str(x) for x in lower_layer)


    return upper_layer_resid, lower_layer_resid

def get_atom_properties (u, dens_type="electron"):
    
    ## assign mass, number of electron or neutron scattering length to each atom types
    
    atom_properties = {}
    
    for name in np.unique (u.atoms.names):
        # Defautl element type in MDAnalysis does not work in some cases
        guessed_element = guess_atom_element(name)
        # for input to pt package make the first letter uppercase and second lowercase
        if len(guessed_element) ==2: guessed_element = guessed_element[0].upper() + guessed_element[1].lower() 
        
        if dens_type == "electron":  
            value = getattr(pt.elements.symbol(guessed_element), "number")
        else:
            value = getattr(pt.elements.symbol(guessed_element), dens_type)
        ## 
        if dens_type == "neutron": value = value.b_c
        
        
        atom_properties [name] = value
        
    return atom_properties
    

def calculate_density (u, group="C114", skip=3, dz=1, begin_frame=0, end_frame=-1, z_lim=40, q0=0.083, N=4, M=4, undulation_correct=1,dens_type="electron"):

     
    # get number of electrons/ mass or neutron sld for each atom type

    atom_properties = get_atom_properties (u, dens_type=dens_type)
    
    shift   =  2*u.dimensions [2] # an arbritrary offset
    
    boxz    =  u.dimensions[2] + shift
     
    Density = np.zeros (( int (boxz/dz), 4))
    
    Density [:,0] = np.linspace (0, boxz,  len(Density))
    
     
     
    Vol_sol = 0
    Vol_lip = 0
    
    # If your water has some other name outside those listed below add that here
    
    possible_water_names = " ".join (x for x in ["H2O", "HOH", "OH2", "HHO", "OHH", "TIP", "T3P", "T4P", "T5P", "SOL", "WAT", "TIP2", "TIP3", "TIP4"])

    # if your ions is not one fo the listed add here
    possible_ion_names   = " ".join (x for x in ["NA", "CL", "K","CA", "SOD","POT","CAL", "CLA", "MG", "NIO", "CXY", "CIO", "LIO", "KIO", "mMg", "nMg"])
    #sol_ion_names =  " ".join (x for x in cg_electron_mass ["SOL"].keys())
  
    membrane = u.select_atoms (f"all and not resname {possible_water_names} {possible_ion_names}")
    z_com = membrane.center_of_mass () [2]

    # get the resids for upper and lower monolayer
    
    upper_layer_resid, lower_layer_resid = get_different_monolayers (u,group) 
    ##################################

    # select group for creating the Undulating reference surfaces for each leaflet based on the provided atom type
    upper_layer_group = u.select_atoms (f"name {group} and resid {upper_layer_resid}")
    lower_layer_group = u.select_atoms (f"name {group} and resid {lower_layer_resid}")
     
    
    ###########################################################################
    mem_resname_string = " ".join (x for x in set(membrane.residues.resnames))    
    Sol_ions = u.select_atoms (f"all and not resname {mem_resname_string}") 
    index=0
    
   ################################################################
        
    for ts in tqdm (u.trajectory [begin_frame:end_frame:skip]):
            
        com_u  =  upper_layer_group.center_of_mass ()
        com_l  =  lower_layer_group.center_of_mass()
        
        # = (z_com_l + z_com_u )/2
    
        com_center = (com_l [2] + com_u [2] )/2
    
        pos_membrane, name_membrane, resname_membrane =  membrane.positions , membrane.names, membrane.resnames #+  shift
         
        pos_membrane [:,2] -= com_center #  center on the center of upper an dlower leaf PO4
        
        # get upper and lower Po4 positions
        R_upper = upper_layer_group.positions
        R_lower = lower_layer_group.positions
        
        # recenter 
        R_upper [:,2] -= com_center #[2]
        R_lower [:,2] -= com_center #[2]
        
        Q, U_filtered = get_fourier_coeff (u.dimensions, R_upper, R_lower, N=N, M=M,q0=q0 ) 
        Cos = []
        #
 
        
        for name, r, resname_ in zip (name_membrane, pos_membrane, resname_membrane):
     
            if undulation_correct:
                
                r_inverse_ft , cos_theta = inverse_transform (r [:2],Q, U_filtered)
                Cos.append (cos_theta)
                #cos_theta=1
                z = cos_theta * (r [2] -r_inverse_ft.real)  + shift

            else: 
                z = r [2] + shift 
                Cos.append(1)
                #print (z)
            #else:  
            #print (cos_theta * (r [2] -r_inverse_ft)-r[2])
            bin_pos = int (z/dz)
           # print (bin_pos)
            if bin_pos > len(Density)-1 or bin_pos<0: continue
    
            Density [bin_pos,2] = Density [bin_pos,2] +  atom_properties [name]  # 10 ppoints

        # Water
        # Water ions
        
        Z_not_membrane, name_not_membrane =  Sol_ions.positions  - com_center   , Sol_ions.names
    
        for name, r in zip (name_not_membrane, Z_not_membrane):

            # one can also apply the undulation correction to the solvent which will mostly affect
            #the solvent near the membrane. But the effect of such correction the overall profiles is
            # negligible. Such a correction slows down the calculation a lot since there are lot of solvent moleucle
            # So we OMIT such corrections
            
            z = r [2] + shift
            bin_pos = int (z/dz)
            
            if bin_pos > len(Density)-1 or bin_pos<0: continue
    
            Density [bin_pos,3] = Density [bin_pos,3] +  atom_properties [name]  # 10 ppoints


        # RNA
 
        Vol_lip = Vol_lip + u.dimensions[1]*u.dimensions[0]*dz/np.average (Cos)
        Vol_sol = Vol_sol + u.dimensions[1]*u.dimensions[0]*dz

    # normalizez
    Density [:, 0] -=  shift  
    Density [:,2] /= Vol_lip
    Density [:,3] /= Vol_sol
    
    ## total density
    Density [:,1] =  Density [:,2] + Density [:,3]
    
    
    # trim the dnesity within the given limits
    if z_lim is None: z_lim = 0.9*u.dimensions [2]/2 # if not set expelicity take 95% of the half the box dimension

    Bool_corrected = (Density [:,0] < z_lim) & (Density [:,0] > -z_lim)
   
    return  Density [Bool_corrected] 


############# input argsmuments

args                = parser.parse_args()
#####################################################

xtc_file            = args.xtc_file
tpr_file            = args.tpr_file
q0                  = args.q0
N                   = args.N
M                   = args.N
skip                = args.skip
dz                  = args.dz
output_filename     = args.output_filename
group               = args.group
z_lim               = args.z_lim
begin_time          = args.begin_time
end_time            = args.end_time
dens_type           = args.dens_type
undulation_correct  = args.undulation_correct
# load trajectory
u = mda.Universe(tpr_file, xtc_file)

begin_frame =   int (begin_time/u.trajectory.dt)

if end_time is None: 
  end_frame = -1
else:
  end_frame   =   int (end_time/u.trajectory.dt)


D = calculate_density (u, group=group,q0=q0, skip=skip, N=N, M=N, dz=dz, z_lim=z_lim, undulation_correct=undulation_correct,dens_type=dens_type,begin_frame=begin_frame, end_frame=end_frame)


#head_group_resid = u.select_atoms (f"name {head_group_atoms}").resids

## write output file
############# Header for the output file
header = "z-distance (Å), total density, lipid density, water+ions density\nThe density units are AMU/Angstrom^3 for mass, e/Angstrom^3 for electron and 10^-6xAngstrom^-2 for neutron sld"

if output_filename is None:
    
    if undulation_correct: 
        np.savetxt (f"{dens_type}_dens_{q0}_{group}_{N}.dat",D, header=header)
    else : 
        np.savetxt (f"{dens_type}_dens_no_undulation_correction.dat", D,header=header)
else:
    
    np.savetxt (output_filename, D, header=header)

