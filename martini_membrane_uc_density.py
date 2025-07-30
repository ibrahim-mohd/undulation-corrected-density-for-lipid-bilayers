import argparse
import json
import pickle
import warnings

import MDAnalysis as mda
import numpy as np
from tqdm import tqdm

# Configure warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser (description="This code calculates transverse mass or electron density profiles for martini membrane or membrane/protein/nucleic acids system. It also perform undulation correction for comparison with expeirments especially form factors")

################ For creating index files of pocket residues ###########################################

parser.add_argument('-f', dest='xtc_file', type=str, default='all.xtc',help='xtc file')
parser.add_argument('-s', dest='tpr_file', type=str, default='em.tpr',help='input tpr file')
parser.add_argument('-j', dest='electron_mass_json', type=str, default='cg_electron.json',help='input json file with CG beads mass or number of electron')
parser.add_argument('-uc', dest='undulation_correct', type=int, default=1, help='if set to 0 no undulation correction is performed')

parser.add_argument('-dz', dest='dz', type=float, default= 1.0, help='slice width in Angstrom')

parser.add_argument('-q0', dest='q0', type=float, default= 0.04, help='filter boundary [nm^-1]')

parser.add_argument('-N', dest='N', type=int, default= 4, help='number of fourier terms')

parser.add_argument('-skip', dest='skip', type=int, default= 1, help='skip frame by this interval')

parser.add_argument('-o', dest='output_filename', type=str, default='./',help='output filename')

parser.add_argument('-group', dest='group', type=str, default='C3A',help='atoms to  consider for fourier coefficient calcualtion')

parser.add_argument('-z_lim', dest='z_lim', type=float, default= None,help='z limit: between -z_lim to +z_lim')
 
        
def get_fourier_coeff(box_dimensions, R_upper, R_lower, N=3, M=3, q0=1.15/10, filter=True):
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

def get_different_monolayers (u, membrane_center_of_mass_z):
    # we obtain the lipids that is int eh top and bottom layers
    # we do that by selecting the head groups/ possible head gorups 
    possible_head_groups = "PO4 NC3 GL1 ROH" # most martini lipids have either PO4 or/and NC3 head groups nad ROH for cholestrol
    upper_layer = np.unique(u.select_atoms (f"name {possible_head_groups} and prop z>{membrane_center_of_mass_z}").residues.resids)
    lower_layer = np.unique(u.select_atoms (f"name {possible_head_groups} and prop z<{membrane_center_of_mass_z}").residues.resids)

    upper_layer_resid = " ".join (str(x) for x in upper_layer)
    lower_layer_resid = " ".join (str(x) for x in lower_layer)

    return upper_layer_resid, lower_layer_resid

    

def get_corrected_density_wth_rna (u, cg_electron_mass, group="C3A", skip=3, dz=1, begin_frame=0, end_frame=-1, z_lim=40, q0=0.083, N=4, M=4, undulation_correct=1):

     
    #n_electrons = {name: electron [name] for name in u.atoms.names}
    
    shift   =  2*u.dimensions [2] # an arbritrary offset
    print (shift)
    boxz    =  u.dimensions[2] + shift
     
    Density = np.zeros (( int (boxz/dz), 4))
    
    Density [:,0] = np.linspace (0, boxz,  len(Density))
    
     
     
    Vol_sol = 0
    Vol_lip = 0

    # select the membrane
 
    sol_ion_names =  " ".join (x for x in cg_electron_mass ["SOL"].keys())
  
    membrane = u.select_atoms (f"all and not resname {sol_ion_names} ION and not name {sol_ion_names}")
    z_com = membrane.center_of_mass () [2]

    # get the resids for upper and lower monolayer
    upper_layer_resid, lower_layer_resid = get_different_monolayers (u, z_com) 
    ##################################

    # select group for creating the Undulating reference surfaces for each leaflet based on the provided atom type
    
    upper_layer_group = u.select_atoms (f"name {group} and resid {upper_layer_resid}")
    lower_layer_group = u.select_atoms (f"name {group} and resid {lower_layer_resid}")
     
    
    ##############
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
    
            Density [bin_pos,2] = Density [bin_pos,2] +  cg_electron_mass [resname_] [name]  # 10 ppoints

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
    
            Density [bin_pos,3] = Density [bin_pos,3] +  cg_electron_mass ["SOL"] [name]  # 10 ppoints


        # RNA
 
        Vol_lip = Vol_lip + u.dimensions[1]*u.dimensions[0]*dz/np.average (Cos)
        Vol_sol = Vol_sol + u.dimensions[1]*u.dimensions[0]*dz

    # normalizez
    Density [:, 0] -=  shift  
    Density [:,2] /= Vol_lip
    Density [:,3] /= Vol_sol

   
     
    
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
undulation_correct  = args.undulation_correct

u = mda.Universe(tpr_file, xtc_file)


if z_lim is None: z_lim = 0.9*u.dimensions [2]/2 # if not set expelicity take 95% of the half the box dimension

#D_cos =  get_density_cos_corrected (u, group=group,q0=q0, skip=20)
with open(args.electron_mass_json, "r", encoding="utf-8") as f: cg_electron_mass = json.load(f)
    
D = get_corrected_density_wth_rna (u, group=group,q0=q0, skip=skip, N=N, M=N, dz=dz, z_lim=z_lim, cg_electron_mass=cg_electron_mass,undulation_correct=undulation_correct)

D [:,1] = D [:,2] + D [:,3] 

############# Header for the output file

header = "z-distance (Å), total density, lipid density, water+ions density\nThe density units are AMU/Angstrom^3 for mass, e/Angstrom^3 for electron and 10^-6xAngstrom^-2 for neutron sld"

 

if undulation_correct: 
    np.savetxt (output_filename.split (".")[0] + f"_{q0}_{group}_{N}.dat",D, header=header)
else : 
    np.savetxt (output_filename.split (".")[0] + "_no_undulation_correction.dat", D,header=header)

