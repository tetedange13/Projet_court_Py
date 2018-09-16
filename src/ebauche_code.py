#!/usr/bin/env python3

import numpy as np
import sys
from Bio.PDB.NACCESS import run_naccess
from Bio.PDB import PDBParser 
import os
import time
from getopt import getopt
import itertools as itr
import multiprocess as mp




def get_acessible_CA(inputFile, path_to_naccess_exe, thresold_ASA):
    """Takes the path to the pdbFile and runs NACCESS on it (ASA calculation)"""
    
    pdb_fileName = os.path.basename(inputFile)
    pdb_id = pdb_fileName.split('.')[0]
    p = PDBParser(QUIET = True)
    structure = p.get_structure(pdb_id, inputFile)
    model = structure[0]
    output_naccess = run_naccess(model, inputFile, 
                                 naccess = path_to_naccess_exe)[0]
    
    #If NACCESS not available:
    #naccessFile = open("../data/6b87.rsa", 'r')
    #output_naccess = naccessFile.readlines()
    #naccessFile.close()
     
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


def get_coord(dict_CA, nb_accessible_CA, nb_tot_CA, inputFile):
    """Get the coordinate of the C-alpha which are accessible to the solvent 
    (i.e. ASA beyond the threesold), from a pdb file (file obj) 
    given as argument. Count also the number of C-alpha"""

    dict_CA_final = dict_CA #To avoid side effects
    
    pdbFile = open(inputFile, 'r')
    pdb_id = os.path.basename(inputFile).split('.')[0]
    
    outFile = open("../results/" + pdb_id + '_out.pdb', 'w')
    
    #To index the resid
    list_resid = [0] * nb_accessible_CA ; i = 0
    #To calculate the center of mass of the protein
    sum_X, sum_Y, sum_Z = 0, 0, 0
    
    resid = 0 #Problems with resid starting at 3, or chains etc in pdbFile
    for line in pdbFile:
        if line.startswith( "ATOM" ):
            outFile.write(line)
        
            if line[12:16].strip() == "CA":
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


    outFile.write("MASTER\n")
    pdbFile.close() ; outFile.close()
    centerOfMass = np.array([sum_X/nb_tot_CA, sum_Y/nb_tot_CA, sum_Z/nb_tot_CA]) 
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



def generate_trigo_arr(precision):
    """Prend un pas (sorte de resolution), decoupe le cercle trigo et renvoie 2
     matrices colonnes des cos et sin de ces angles"""
     
    #Prop of trigo functions trigo used:
    #cos(pi-x) = -cos(x) ; sin(pi-x) = sin(x) <=> N-O dial
    #cos(pi+x) = -cos(x) ; sin(pi+x) = -sin(x) <=> S-O dial
    #cos(-x) = cos(x) ; sin(-x) = -sin(x) <=> S-E dial
    
    size_arr = precision+1
    arr_cos = np.zeros( size_arr, dtype=float )
    arr_sin = np.zeros( size_arr, dtype=float )
    
    
    #We start with filling with evident values (not zero) of cos et sin:
    #En 0:
    arr_cos[0] = 1.0
    #En pi:
    arr_cos[precision] = -1.0
    
    if precision%2 == 0: #If the precision is even, there is also pi/2
        idx_pi_over_2 = int( size_arr/2 )
        #En pi/2:
        arr_sin[idx_pi_over_2] = 1.0 
            
            
    #N-E DIAL (will serve as a reference to fill the other dial):
    i=1
    for angle in np.arange(np.pi/precision, 
                           np.pi/2 - 1/precision, 
                           np.pi/precision):
    #J'ai trouve le "-1/pas" un peu empiriquement => A VOIR...
        arr_cos[i] = np.cos(angle)
        arr_sin[i] = np.sin(angle)
        i += 1     
    
    #N-O dial:  
    nb_to_calc = i - 1 #Number of values to calculate
    
    #Depending on the parity of the precision, the start changes:
    if precision%2 == 0:  
        start = idx_pi_over_2 + 1
    else:
        start = i 
        
    end = start + nb_to_calc
    
    arr_cos[start:end] = -arr_cos[1:i][::-1]
    arr_sin[start:end] = arr_sin[1:i][::-1]
    
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


#To simplify the syntax, in the next functions, the use of "CA" in any variable
#will refer to "all the CA which are enough accessible to the solvent"

    
def start_position_plane(unit_vect, dict_CA, nb_CA, list_resid):
    """Takes an array containing the coord of a unit vect.
    1) Puts the plane orthogonal to the unit vect given as argument, far enough 
    from the protein (to be completely outside of it)
    2) Calculates the distances between all the CA of the protein and this plane
    ( values are stocked in a array of size (nb_accessible_CA, 1) )
    3) Determine both the closer and the further CA from the given plane, in 
    order to calculate the number of steps (=nb_slides) needed to cross the whole
    protein along the given axis.
    4) Translate the distance values, in order to position the plan just under
    the closer CA
    4) Returns this nb_slides and the (transformed) array of distances, which 
    will be used to slide along the axis (function next next)"""

    #Start plane, far from the protein:
    initial_dist = 1000
    initial_plane = np.hstack( (initial_dist*unit_vect, -initial_dist**2) )

    #Will contain dist between each point and the plane:
    start_dist_points_arr = np.zeros((nb_CA, 1)) 
    
    i = 0 
    for resid in dict_CA.keys():
        start_dist_points_arr[i, ] = dist_point_plane( dict_CA[resid], \
                                                initial_plane )
        list_resid[i] = resid
        i += 1

    start_dist_plane = np.floor( initial_dist - np.min(start_dist_points_arr) )
    end_dist_plane = np.floor( initial_dist - np.max(start_dist_points_arr) ) 
    #=> Que faire si la dist entre le plus pres et le plus loi es trop petite ??
    
    #We transform the dist array, to adapt to the well positionned plane:
    start_dist_points_arr -= (initial_dist - start_dist_plane)
    nb_slides = int( start_dist_plane - end_dist_plane + 15 ) #OUI je pense

    
    #For further improval, we need the start dist for the plan (last arg):
    return ( nb_slides, start_dist_points_arr, start_dist_plane )



def freq_hydrophob(start_dist_points_arr, nb_CA, dict_CA, list_resid):
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
        if -15 <= start_dist_points_arr[i, ] <= 0: #For the aa in the slice
            nb_accessible += 1
            
            resid = list_resid[i]
            aa = dict_CA[resid]["resName"]
            
            if aa in list_hydrophob_aa:
                nb_hydrophob_accessible += 1
    
    #if nb_hydrophob_accessible/nb_accessible > 0.3:
        #print("NB_HYDROPHB= ",nb_hydrophob_accessible, "NB_ACCESS= ",nb_accessible)
        #print("FREQ = ", nb_hydrophob_accessible/nb_accessible)
    return nb_hydrophob_accessible/nb_accessible



def sliding_slice(nb_slides, start_dist_points_arr, nb_CA, dict_CA, list_resid):
    """Takes a slice positionned at its start point, by taking 
    its start_dist_points_arr. 
    Slides this slice nb_slides times, by decrementing the dist_arr of 1A.
    At each step, the frequency of hydrophobic aa among the accessible aa is
    calculated, in order to have the mean of the frequency of hydrophobic aa
    along the current direction
    
    REMARK: The different values of frequency along a given direct are not kept
    individually, in order to reduce the cost in memory.
    That is why, after found the good normal vect to the membrane, we will have 
    to slide one last time along this direction, in order to well position the
    membrane"""

    dist_arr = start_dist_points_arr
    sum_freq_hydrophob = 0
    step_slides = 1
    
    for r in range(0, nb_slides, step_slides):
        sum_freq_hydrophob += freq_hydrophob( dist_arr, 
                                              nb_CA, 
                                              dict_CA, 
                                              list_resid )
        dist_arr -= step_slides #We slide the slice of 1A
        

    return sum_freq_hydrophob/nb_slides
         


def dist_point_origin(point):
    """Calculates the distance between a point (given by its coordinates inside
    a dict) and the origin of the system"""
    
    arr_point = dict_to_arr(point)
    
    return np.sqrt( arr_point @ arr_point.T )



def improve_mb_position(best_direction, dict_CA, nb_CA, 
                        list_resid, best_start_dist_plane):
    #Je l'ai un peu optimisée, pour éviter de stocker en memoire une matrice de
    #freq, que je dois ensuite parcourir pour trouver le max avec np.argmax()
    #Par contre, le repositionnement du plan au debut, il est pas tres 
    #optimal, puisque que tu dois recalculer la dist entre le plan et tous les 
    #points...

    """
    1) Slides along the given direction, in order to find the position which has
    the best value of freq_hydrophob
    2) Makes the membrane grows on each side (0.1 by 0.1A), while this increases
    the freq_hydrophob inside the slice, in order to really to find the optimal
    thickness of the membrane"""
    
    nb_slides, start_dist_points_arr = start_position_plane( best_direction, 
                                                             dict_CA, 
                                                             nb_CA, 
                                                             list_resid ) [0:2]
    
    dist_arr = start_dist_points_arr
    
    current_max_freq = freq_hydrophob(dist_arr, nb_CA, dict_CA, list_resid)
    dist_arr -= 1 ; idx_max_freq = 0

    for r in range(1, nb_slides):
        new_freq = freq_hydrophob(dist_arr, nb_CA, dict_CA, list_resid)
      
        if new_freq > current_max_freq:
            current_max_freq = new_freq
            idx_max_freq = r
            
        dist_arr -= 1 #We slide the slice of 1A
    
        
    return best_start_dist_plane - idx_max_freq


def draw_point(point, num, outFile):
    well_formated = "{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}\n"

    outFile.write( well_formated.format( "HETATM", num, "N", "DUM", num,
                                         point[0],
                                         point[1],
                                         point[2] ))
    

def draw_axis(start, end, step, best_direction, centerOfMass, outFile):
    """Draws an axis from a direction, given by a unit vector"""               
    for r in range(start, end, step):
        point_to_write = r * best_direction + centerOfMass
        draw_point(point_to_write, r, outFile)

        
                                        
def draw_plane(dist_best_plane, idx_best_angles, precision, 
               centerOfMass, outFile):
    gap = 5 #Distance between two consecutive points
    limit_angle = 3*np.pi/7
    nb_points = 5
    unit_angle = limit_angle/nb_points
    
    idx_theta, idx_phi = idx_best_angles
    theta, phi = idx_best_angles * np.pi / precision
    
    for i in range(1, nb_points):
        angle = np.arctan(i * gap / dist_best_plane) 
        H = dist_best_plane / np.cos(phi + angle)
        
        for j in range(8):
            X_point = H * np.sin(phi + angle) * np.cos(theta + j*np.pi/4)
            Y_point = H * np.sin(phi + angle) * np.sin(theta + j*np.pi/4)
            Z_point = H * np.cos(phi + angle)
            point = np.array( [X_point, Y_point, Z_point] ) + centerOfMass
            draw_point(point, i, outFile)

            
            
def draw_planes(dist_best_plane, trigo_phi, trigo_theta, centerOfMass,
                arr_trigo_phiShift, arr_trigo_thetaShift, 
                outFile, mode="hollow"):
    #cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
    #sin(a + b) = sin(a)*cos(b) + sin(b)*cos(a)
    
    cos_phi, sin_phi = trigo_phi ; cos_theta, sin_theta = trigo_theta
    arr_cos_phiShift, arr_sin_phiShift = arr_trigo_phiShift
    arr_cos_thetaShift, arr_sin_thetaShift = arr_trigo_thetaShift

    if mode == "plain":
        width_circle = 3 #Number of points in the radium of the external circle
        start = nb_points - width_circle
    elif mode == "hollow":
        start = 1
    
    for i in range(start, nb_points):
        cos_phiShift, sin_phiShift = arr_cos_phiShift[i], arr_sin_phiShift[i]
        
        cos_phiShifted = cos_phi*cos_phiShift - sin_phi*sin_phiShift
        sin_phiShifted = sin_phi*cos_phiShift + cos_phi*sin_phiShift
        
        H = dist_best_plane / cos_phiShifted     
        period = len(arr_cos_phiShift)
        
        for j in range(2*period):
            cos_thetaShift = arr_cos_thetaShift[j]
            sin_thetaShift = arr_sin_thetaShift[j]
            cos_thetaShifted = cos_theta*cos_thetaShift-sin_theta*sin_thetaShift
            sin_thetaShifted = sin_theta*cos_thetaShift+cos_theta*sin_thetaShift
        
            X_point = H * sin_phi_shifted * cos_thetaShifted
            Y_point = H * sin_phi_shifted * sin_thetaShifted
            Z_point = H * cos_phi_shifted
            point = np.array( [X_point, Y_point, Z_point] ) + centerOfMass
            draw_point(point, i, outFile)
            


def calc_coord_plane(inputFile, best_direction, dist_best_plane, centerOfMass,
                     arr_cos, arr_sin, mode="hollow"):
    """To be able to have a nice PyMol output of our plan, we need to determine
    the coordinates of the 4 corners of our plane 
    (represented as a 3D rectangle)"""
    
    pdb_id = os.path.basename(inputFile).split('.')[0]
    outFile = open("../results/" + pdb_id + "_out.pdb", 'a')
    
    
    #Ici, on va utiliser la formule suivante:
    # cos(phi + pi/4) = R / H
    # <==> H = R / cos(phi + pi/4)
    
    #x = r * sin(phi) * cos(theta)
    #y = r * sin(phi) * sin(theta)
    #z = r * cos(phi)
    
    #Draw the axis orthogonal to the membrane
    draw_axis(0, dist_best_plane, 5, best_direction, centerOfMass, outFile)
    draw_axis(-25, 25, 5, best_direction, centerOfMass, outFile)
        
    idx_theta, idx_phi = idx_best_angles
    trigo_phi = np.array([ arr_cos[idx_phi], arr_sin[idx_phi] ])
    trigo_theta = np.array([ arr_cos[idx_theta], arr_sin[idx_theta] ])
    
    limit_angle =  2*np.pi/7
    nb_points = 10
    phiShift = limit_angle/nb_points #Angle between two consecutive points
    #arr_trigo_phiShift = generate_trigo_arr(phiShift)
    
    density = 8 ; thetaShift = np.pi/density
    #arr_trigo_thetaShift = generate_trigo_arr(theta_Shift)

    #draw_plane(dist_best_plane, idx_best_angles, precision, centerOfMass, outFile)
    #draw_planes(dist_best_plane+15, trigo_phi, trigo_theta, centerOfMass,
                #arr_trigo_phiShift, arr_trigo_thetaShift, outFile)   
    
    
    #We write a pymol_file, which can be executed automaticly:
    pymol_file = open("cmd_pymol.pml", 'w')
    
    loading = "cmd.load('" + str(outFile.name) + "')\n"
    bg = "cmd.bg_color('white')\n"
    sticks = "cmd.show('sticks')\n"
    spheres = "cmd.show('spheres', 'HETATM')\n"
    
    for ligne in (loading, bg, sticks, spheres):
        pymol_file.write(ligne)
    
    pymol_file.close()
    
    outFile.close()
    
    
    
    
    #HETATM 2548  N   DUM  2548     -16.000  -6.000 -15.700
    #"{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(...)
    
    return

 
def parallelized_fun(tupl):
    i,j = tupl
    unit_vect = arr_unit_vect[:, i, j]
            
    #4-1) We position the plan well:
    tuple_returned = start_position_plane( unit_vect, dict_CA, 
                                               nb_accessible_CA, list_resid )
    nb_slides, start_dist_points_arr, start_dist_plane = tuple_returned
    start_dist_plane_arr[i, j] = start_dist_plane
    arr_nb_slides[i, j] = nb_slides
        
    #4-2) We slide along the current direction:   
    return sliding_slice( nb_slides, start_dist_points_arr, 
                                             nb_accessible_CA, dict_CA, 
                                             list_resid )      
   


#MAIN

start_time = time.time()

#We get the arguments given:
path_to_naccess_exe = 'naccess'
opts, args = getopt(sys.argv[1:], "hi:p:", ["naccess=", "precision="])

for opt, arg in opts:
    if opt == "-h":
        print("Ceci est un message d'aide")
        sys.exit()
        
    elif opt == "-i":
        inputFile = arg    
        
    elif opt == "--naccess":
        path_to_naccess_exe = arg
    
    elif opt in ("-p", "--precision"):
        precision = int(arg)       


# 1) We extract the needed data:
thresold_ASA = 30 #Considered accessible if ASA up to thresold

dict_CA, nb_accessible_CA, nb_tot_CA = get_acessible_CA( inputFile, 
                                                         path_to_naccess_exe, 
                                                         thresold_ASA)


dict_CA, list_resid, centerOfMass = get_coord( dict_CA, nb_accessible_CA, 
                                               nb_tot_CA, inputFile)


# 2) We transform the coordinates to set the center of mass as the origin:
dict_CA = transform_coord(dict_CA, centerOfMass)


# 3) We generate the  arrays:
#Formulas, from polar to cartesian:
#x = r * sin(phi) * cos(theta)
#y = r * sin(phi) * sin(theta)
#z = r * cos(phi)

#These arrays need to be created only one:
arr_cos_theta, arr_sin_theta = generate_trigo_arr(precision)
arr_cos_phi, arr_sin_phi = arr_cos_theta, arr_sin_theta

#Use of matricial product to calculate each of the x, y and z matrixes
#of the unit vector dividing half of the the 3D space
#(necessity to reshape the different matrixes, to make the product)
size_arr = precision + 1
vectArr_X = arr_cos_theta[:, np.newaxis] @ arr_sin_phi[np.newaxis, :]
vectArr_Y = arr_sin_theta[:, np.newaxis] @ arr_sin_phi[np.newaxis, :]
vectArr_Z = np.tile( arr_cos_phi, (size_arr, 1) )

#REMARK: The vectArr_Z has been created by repeating the line of cos(phi) as many
#times there are values for the theta angle

#We merge the 3 arrays in a single one, to simplify its passing as an argument: 
arr_unit_vect = np.zeros( (3, size_arr, size_arr), dtype = float )
arr_unit_vect[0, :, :] = vectArr_X
arr_unit_vect[1, :, :] = vectArr_Y
arr_unit_vect[2, :, :] = vectArr_Z


# 4) Loop To calculate the mean freq_hydrophob on all the directions:
arr_mean_freq = np.zeros( (size_arr, size_arr) )
start_dist_plane_arr = np.zeros( (size_arr, size_arr) , dtype = int)
arr_nb_slides = np.zeros( (size_arr, size_arr) , dtype = int)

#for i in range(size_arr):
#    for j in range(size_arr):


#Parallelization
input = ((i, j) for i, j in itr.combinations_with_replacement(range(size_arr), 2))
p = mp.Pool()
flat_mean_freq = np.array(p.map(parallelized_fun, input))
arr_mean_freq = flat_mean_freq.reshape( (size_arr, size_arr) )
p.close()
p.join()
print(results)
        

print("MAX_FREQ = ", np.max(arr_mean_freq))
print("\nESPACE\n")
print(np.nonzero(arr_mean_freq - np.max(arr_mean_freq)) == 0)

# 5) Determination of the indexes of this maximum value 
#inside the arr_mean_freq:
idx_max_flat = np.argmax(arr_mean_freq)
idx_max_goodShape = np.unravel_index( idx_max_flat, (size_arr, size_arr) )
idx_theta = idx_max_goodShape[0] ; idx_phi = idx_max_goodShape[1]
idx_best_angles = np.array( [idx_theta, idx_phi], dtype=int )

best_direction = arr_unit_vect[:, idx_theta, idx_phi]
theta = str(idx_theta) + 'pi/' + str(precision)
phi = str(idx_phi) + 'pi/' + str(precision)
print("THETA = ", theta, " ; ", "PHI = ", phi)


best_start_dist_plane = start_dist_plane_arr[idx_theta, idx_phi]

ma_direction = arr_unit_vect[:, 3, 4]
print('\n ESPACE \n')
dist_best_plane = improve_mb_position( ma_direction, dict_CA, nb_accessible_CA, 
                                       list_resid, best_start_dist_plane )


# 6) Manage the output:

calc_coord_plane(inputFile, best_direction, dist_best_plane, centerOfMass,
                     arr_cos_phi, arr_sin_phi)

print("RUNTIME = ", time.time() - start_time)












