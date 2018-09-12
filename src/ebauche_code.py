#!/usr/bin/env python3

import numpy as np
import sys
from Bio.PDB.NACCESS import run_naccess
from Bio.PDB import PDBParser 
import os



def get_acessible_CA(pdb_id, path_to_naccess_exe, thresold_ASA):
    """Takes the pdb id (as a str), runs NACCESS on it (ASA calculation)"""
    
    pdb_fileName = "../data/" + pdb_id + ".pdb"
    p = PDBParser()
    structure = p.get_structure(pdb_id, pdb_fileName)
    model = structure[0]
    output_naccess = run_naccess(model, pdb_fileName, 
                                 naccess = path_to_naccess_exe)[0]
    
    dict_CA = {} #Contains all useful info for CA
    nb_accessible_CA = 0
    resid = 0
    
    for line in output_naccess:
        splitted_line = line.split()
        type_line = splitted_line[0]
        resName = splitted_line[1]
        
        if type_line == "RES" and len(resName) == 3:
            resid += 1
            numRes = int( splitted_line[3] )
            relASA = float(splitted_line[5])
            chain = splitted_line[2]
            
            if relASA >= thresold_ASA:
                nb_accessible_CA += 1
                dict_CA[resid] = { 'resName': resName,
                                   'chain': chain,
                                   'numRes': numRes, 
                                   'relASA': relASA }
    
    
    #The resid correspond here to the total number of CA             
    return ( dict_CA, nb_accessible_CA, resid )


def get_coord(dict_CA, nb_accessible_CA, nb_tot_CA, pdbFile):
    """Get the coordinate of the C-alpha which are accessible to the solvent 
    (i.e. ASA beyond the threesold), from a pdb file (file obj) 
    given as argument. Count also the number of C-alpha"""
    
    dict_CA_final = dict_CA #To avoid side effects
    
    #To index the resid
    list_resid = [0] * nb_accessible_CA ; i = 0
    #To calculate the center of mass of the protein
    sum_X, sum_Y, sum_Z = 0, 0, 0
    
    resid = 0 #Problems with resid starting at 3, or chains etc in pdbFile
    for line in pdbFile:
        if line.startswith( "ATOM" ) and line[12:16].strip() == "CA":
            resid += 1
            X_coord = float( line[30:38].strip() )
            Y_coord = float( line[38:46].strip() )
            Z_coord = float( line[46:54].strip() )
            
            #To calculate the center of mass of the protein:
            sum_X += X_coord ; sum_Y += Y_coord ; sum_Z += Z_coord

            if resid in dict_CA.keys(): #If considered as accessible
                list_resid[i] = resid ; i += 1
                
                dict_CA_final[resid]['x'] = X_coord 
                dict_CA_final[resid]['y'] = Y_coord 
                dict_CA_final[resid]['z'] = Z_coord  
    

    centerOfMass = [ sum_X/nb_tot_CA, sum_Y/nb_tot_CA, sum_Z/nb_tot_CA ] 
    return (dict_CA_final, list_resid, centerOfMass)


 
def transform_coord(dict_coord, centerOfMass):
    """Takes the coordinates of the CA of the protein, its center of mass and
    returns a dict with transformed coord, i.e. with the center of mass set as 
    the origin of the system."""
 
    transformed_dict = dict_coord #To avoid side effects
    for resid in dict_coord.keys():
        transformed_dict[resid]['x'] -= centerOfMass[0]
        transformed_dict[resid]['y'] -= centerOfMass[1]
        transformed_dict[resid]['z'] -= centerOfMass[2]
        
        
    return transformed_dict



def generate_sinCos_arr(precision):
    """Prend un pas (sorte de resolution), decoupe le cercle trigo et renvoie 2
     matrices colonnes des cos et sin de ces angles"""
     
    #Prop of trigo functions trigo used:
    #cos(pi-x) = -cos(x) ; sin(pi-x) = sin(x) <=> N-O dial
    #cos(pi+x) = -cos(x) ; sin(pi+x) = -sin(x) <=> S-O dial
    #cos(-x) = cos(x) ; sin(-x) = -sin(x) <=> S-E dial
    
    
    #We 1st disinguish the cases where we have to deal with theta (from 0 to 2pi)
    #Or phi (from à to pi):
   
    arr_cos = np.zeros((precision+1, 1), dtype=float)
    arr_sin = np.zeros((precision+1, 1), dtype=float)
    

    #We start with filling with evident values (not zero) of cos et sin:
    #En 0:
    arr_cos[0, 0] = 1.0
    #En pi:
    arr_cos[precision, 0] = -1.0
    
    if precision%2 == 0: #If the precision is even, there is also pi/2
        idx_pi_over_2 = int( (precision+1)/2 )
        #En pi/2:
        arr_sin[idx_pi_over_2, 0 ] = 1.0 
            
            
    #N-E DIAL (will serve as a reference to fill the other dial):
    i=1
    for angle in np.arange(np.pi/precision, 
                           np.pi/2 - 1/precision, 
                           np.pi/precision):
    #J'ai trouve le "-1/pas" un peu empiriquement => A VOIR...
        arr_cos[i, 0] = np.cos(angle)
        arr_sin[i, 0] = np.sin(angle)
        i += 1     
    
    #N-O dial:  
    nb_to_calc = i - 1 #Number of values to calculate
    
    #Depending on the parity of the precision, the start changes:
    if precision%2 == 0:  
        start = idx_pi_over_2 + 1
    else:
        start = i 
        
    end = start + nb_to_calc
    
    arr_cos[start:end, 0] = -arr_cos[1:i, 0][::-1]
    arr_sin[start:end, 0] = arr_sin[1:i, 0][::-1]
    
    return (arr_cos, arr_sin)        
    


def dict_to_arr(point):
    """Takes a point as a dict and convert it in a array(3,1)"""
    
    return np.array( [ point['x'], point['y'], point['z'] ] )
        


def dist_point_plane(point, plane):
    """Prend un point provenant du pdb (donc un dict)
    et un plan = (a, b; c; d), avec:
        -- (a; b; c)  = (x; y; z) du point au bout d'un vecteur 
        -- et d = -Rcarre"""
    
    arr_point = np.hstack( (dict_to_arr(point), 1) )
        
    numer = abs( arr_point @ plane.T )
    denom = np.sqrt( plane[0:3, ] @ plane[0:3, ].T )
    return numer/denom


#To simplify the syntax, in the next functions, the use of "nb_CA" will refer to
# "the number of CA which are enough accessible to the solvent"

    
def start_position_plane(unit_vect, coord_pdbFile, nb_CA, list_resid):
    """Takes an array containing the coord of a unit vect.
    1) Puts the plane orthogonal to the unit vect given as argument, far enough 
    from the protein (to be completely outside of it)
    2) Calculates the distances between all the CA of the protein and this plane
    ( values are stocked in a array of size (nb_accessible_CA, 1) )
    3) Determine both the closer and the further CA from the given plane, in 
    order to calculate the number of steps (=nb_slides) needed to cross the whole
    protein along the given axis.
    4) Translate the distance values, in order to positionne the plan just under
    the closer CA
    4) Returns this nb_slides and the (transformed) array of distances, which 
    will be used to slide along the axis (function next next)"""

    #Start plane, far from the protein:
    initial_dist = 1000
    initial_plane = np.hstack( (initial_dist*unit_vect, -initial_dist**2) )
    #Will contain dist between each point and the plane:
    start_dist_arr = np.zeros((nb_CA, 1)) 
    
    i = 0 
    for resid in coord_pdbFile.keys():
        start_dist_arr[i, ] = dist_point_plane( coord_pdbFile[resid], \
                                                initial_plane )
        list_resid[i] = resid
        i += 1

    dist_closer = np.floor( initial_dist - np.min(start_dist_arr) )
    dist_further = np.floor( initial_dist - np.max(start_dist_arr) ) 
    #=> Que faire si la dist entre le plus pres et le plus loi es trop petite ??
    
    #We transform the dist array, to adapt to the well positionned plane:
    start_dist_arr -= (initial_dist - dist_closer)
    nb_slides = int( dist_closer - dist_further + 15 ) #OUI je pense
    
    return ( nb_slides, start_dist_arr )



def freq_hydrophob(start_dist_arr, nb_CA, dict_CA, list_resid):
    """Read the dist_arr, in order to find which one are in the current slice.
    Then goes in the NACCESS output file, to know if it is an hydrophobic
    or not
    Return the frequence of hydrophobic among all accessible aa in the current
    slice"""

    nb_accessible, nb_hydrophob_accessible = 0, 0
    #The  aa considered as hydrophobic in the article:
    list_hydrophob_aa = ["PHE", "GLY", "ILE", "LEU", 
                         "MET", "VAL", "TRP", "TYR"] 
    
    for i in range(nb_CA):
        if -15 <= start_dist_arr[i, ] <= 0: #For the aa in the slice
            nb_accessible += 1
            
            resid = list_resid[i]
            aa = dict_CA[resid]["resName"]
            
            if aa in list_hydrophob_aa:
                nb_hydrophob_accessible += 1
            
    return nb_hydrophob_accessible/nb_accessible



def sliding_slice(nb_slides, start_dist_arr, nb_CA, dict_CA, list_resid):
    """Takes a slice positionned at its start point, by taking 
    its start_dist_arr. 
    Slides this slice nb_slides times, by decrementing the dist_arr of 1A.
    At each step, the frequency of hydrophobic aa among the accessible aa is
    calculated, in order to have the mean of the frequency of hydrophobic aa
    along the current direction"""

    dist_arr = start_dist_arr
    sum_freq_hydrophob = 0
    for r in range(nb_slides):
        print( freq_hydrophob(dist_arr, nb_CA, dict_CA, list_resid) )
        sum_freq_hydrophob += freq_hydrophob(dist_arr, 
                                             nb_CA, 
                                             dict_CA, 
                                             list_resid)
        dist_arr -= 1 #We slide the slice of 1A
        
    return sum_freq_hydrophob/nb_slides
    


def is_in_slice(point, planeDown, planeUp):
    """Determine if a given C-alpha (aka point with (x;y;z) coord into a dict) 
    is inside a given slice of 1A wide (between planeI and planeI_plus1)"""
    
    if dist_point_plane(point, planeDown) < 1 \
                                and dist_point_plane(point, planeUp) < 1:
        return 1    
    
    else:
        return 0                                
    
    
   
def fonction_pple (r, point):
    size_arr_theta = np.shape(arr_cos_theta)[0]  
    vectArr_X_down, vectArr_X_up = r * arr_cos_theta @ arr_sin_phi.T, \
                           (r+1) * arr_cos_theta @ arr_sin_phi.T
    vectArr_Y_down, vectArr_Y_up = r * arr_sin_theta @ arr_sin_phi.T, \
                           (r+1) * arr_sin_theta @ arr_sin_phi.T
    vectArr_Z_down, vectArr_Z_up = r * np.tile( arr_cos_phi.T, (size_arr_theta, 1) ), \
                           (r+1) * np.tile( arr_cos_phi.T, (size_arr_theta, 1) )  
    
    taille_i, taille_j = np.shape(vectArr_X_down)
    for i in range(precision + 1):
        for j in range(precision + 1):
            planeDown = np.array( [vectArr_X_down[i,j],
                                   vectArr_Y_down[i,j],
                                   vectArr_Z_down[i,j],
                                   -r**2 ] )
            planeUp = np.array( [vectArr_X_up[i,j],
                                   vectArr_Y_up[i,j],
                                   vectArr_Z_up[i,j],
                                   -(r+1)**2 ] )
            print( is_in_slice(point, planeDown, planeUp) )
         



#MAIN

arg_cmd = sys.argv


# 1) We extract the needed data:
pdb_id = "6b87"
thresold_ASA = 20 #Considered accessible if ASA up to thresold

if len(arg_cmd) == 2:
    path_to_naccess_exe = 'naccess'
else:
    path_to_naccess_exe = arg_cmd[2]

dict_CA, nb_accessible_CA, nb_tot_CA = get_acessible_CA( pdb_id, 
                                                         path_to_naccess_exe, 
                                                         thresold_ASA)
print("NB_CA_ACCESSIBLE = ", nb_accessible_CA)
print("NB_CA_TOT = ", nb_tot_CA)

pdbFile = open('../data/' + pdb_id + ".pdb", 'r')
dict_CA, list_resid, centerOfMass = get_coord( dict_CA, 
                                               nb_accessible_CA, 
                                               nb_tot_CA,
                                               pdbFile)
pdbFile.close()


# 2) We transform the coordinates to set the center of mass as the origin:
dict_CA = transform_coord(dict_CA, centerOfMass)


# 3) We generate the  arrays corresponding respecsinCos_arr
#Formulas, from polar to cartesian:
#x = r * sin(phi) * cos(theta)
#y = r * sin(phi) * sin(theta)
#z = r * cos(phi)

precision = int(arg_cmd[1])


#These arrays need to be created only one:
arr_cos_theta, arr_sin_theta = generate_sinCos_arr(precision)
arr_cos_phi, arr_sin_phi = arr_cos_theta, arr_sin_theta

#Use of matricial product to calculate each of the x, y and z matrixes
#of the unit vector dividing the 3D space
#(necessity to reshape the phi matrixes, to make the product)
vectArr_X = arr_cos_theta @ arr_sin_phi.T
vectArr_Y = arr_sin_theta @ arr_sin_phi.T
vectArr_Z = np.tile( arr_cos_phi.T, ((precision+1), 1) )

#REMARK: The vectArr_Z has been created by repeating the line of cos(phi) as many
#times there are values for the theta angle


#TOUT CA, CA DOIT ETRE BOUCLEH, SUR L'ENSEMBLE DES DIRECTIONS:

mon_vect = np.array( [vectArr_X[4, 4], vectArr_Y[4, 4], vectArr_Z[4, 4]] )
nb_slides, start_dist_arr = start_position_plane(mon_vect, 
                                                 dict_CA, 
                                                 nb_accessible_CA,
                                                 list_resid)



mean_this_direction = sliding_slice(nb_slides, 
                       start_dist_arr, 
                       nb_accessible_CA, 
                       dict_CA, 
                       list_resid)
print(mean_this_direction)













