#!/usr/bin/env python3

from src.parse_files import *
from src.geometry import *
import sys
import time
from getopt import getopt
import itertools as itr
import multiprocess as mp




#To simplify the syntax, in the next functions, the use of "CA" in any variable
#will refer to "all the CA which are enough accessible to the solvent"

    
def start_position_plane(unit_vect, dict_CA, nb_CA, list_resid):
    """Takes an array containing the coord of a unit vect.
    1) Puts the plane orthogonal to the unit vect given as argument, far enough 
    from the protein (to be completely outside of it)
    2) Calculates the distances between all the CA of the protein and this plane
    ( values are stocked in a array of size (nb_accessible_CA, 1) )
    3) Determine both the closer and the further CA from the given plane, in 
    order to calculate the number of steps (=nb_slices) needed to cross the whole
    protein along the given axis.
    4) Translate the distance values, in order to position the plan just under
    the closer CA
    4) Returns this nb_slices and the (transformed) array of distances, which 
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

    start_dist_plane = np.floor(initial_dist - np.min(start_dist_points_arr)) +1
    end_dist_plane = np.floor( initial_dist - np.max(start_dist_points_arr) ) 
    
    #We transform the dist array, to adapt to the well positionned plane:
    start_dist_points_arr -= (initial_dist - start_dist_plane)
    nb_slices = int( start_dist_plane - end_dist_plane ) - 15

    
    #For further improval, we need the start dist for the plan (last arg):
    return ( nb_slices, start_dist_points_arr, start_dist_plane )



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
        if 0 <= start_dist_points_arr[i, ] <= 15: #For the aa in the slice
            nb_accessible += 1
            
            resid = list_resid[i]
            aa = dict_CA[resid]["resName"]
            
            if aa in list_hydrophob_aa:
                nb_hydrophob_accessible += 1
    
    
    return nb_hydrophob_accessible/nb_accessible



def sliding_slice(nb_slices, start_dist_points_arr, nb_CA, dict_CA, list_resid):
    """Takes a slice positionned at its start point, by taking 
    its start_dist_points_arr. 
    Slides this slice nb_slices times, by decrementing the dist_arr of 1A.
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
    step_slides = 5 #Seems to be the optimal value for this parameter
    
    for r in range(0, nb_slices, step_slides):
        sum_freq_hydrophob += freq_hydrophob( dist_arr, nb_CA, 
                                              dict_CA, list_resid )       
        dist_arr -= step_slides #We slide the slice of 1A
    
    
    #Needed to add the square of nb_slices to avoid the biais on nb_slices
    return nb_slices**2 * sum_freq_hydrophob/nb_CA
         


def parallelized_fun(tupl):
    """To calculate the mean freq_hydrophob on all the directions (in a 
    parallelized way)"""
    
    i,j = tupl
    unit_vect = arr_unit_vect[:, i, j]
    
            
    #4-1) We position the plan well:
    tuple_returned = start_position_plane( unit_vect, dict_CA, nb_accessible_CA,
                                           list_resid )
    nb_slices, start_dist_points_arr, start_dist_plane = tuple_returned
        
    #4-2) We slide along the current direction:   
    return ( start_dist_plane, sliding_slice(nb_slices, start_dist_points_arr, 
                                             nb_accessible_CA, dict_CA,
                                             list_resid) )         



def improve_mb_position(best_direction, dict_CA, nb_CA, list_resid):
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
    
    nb_slices, start_dist_points_arr = start_position_plane( best_direction, 
                                                             dict_CA, 
                                                             nb_CA, 
                                                             list_resid ) [0:2]
    
    dist_arr = start_dist_points_arr
    
    current_max_freq = freq_hydrophob(dist_arr, nb_CA, dict_CA, list_resid)
    dist_arr -= 1 ; idx_max_freq = 0

    for r in range(1, nb_slices):
        new_freq = freq_hydrophob(dist_arr, nb_CA, dict_CA, list_resid)

        if new_freq > current_max_freq:
            current_max_freq = new_freq
            idx_max_freq = r
            
        dist_arr -= 1 #We slide the slice of 1A
    

    return idx_max_freq


def draw_point(point, num, outFile):
    well_formated = "{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}\n"

    outFile.write( well_formated.format( "HETATM", num, "N", "DUM", num,
                                         point[0],
                                         point[1],
                                         point[2] ))
    

def draw_axis(start, end, best_direction, centerOfMass, outFile):
    """Draws an axis from a direction, given by a unit vector""" 
    
    step = 5
    if start == 0 and end < 0: step *= -1    
                  
    for r in range(start, end, step):
        point_to_write = r * best_direction
        draw_point(point_to_write, r, outFile)

        
                                        
def draw_plane_2(dist_best_plane, idx_best_angles, precision, 
               centerOfMass, outFile):
    gap = 5 #Distance between two consecutive points
    limit_angle = 3*np.pi/7
    nb_points = 5 * dist_best_plane
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



def draw_plane(best_direction, dist_best_plane, centerOfMass, outFile):
    """Calculate the coordinates of the points belonging to the planes which
    are normal to the best direction and draw these planes by writing (appending)
    DUM atoms inside a pdb file"""
    
    best_vect = best_direction * dist_best_plane
    
    produit = np.dot( np.array([0.0, 0.0, 1.0])[np.newaxis, :],
                      best_vect[:, np.newaxis] )[0][0]
    cos_transform = produit / dist_best_plane                  
    sin_transform = np.sqrt(1 - cos_transform**2)
    transform_arr = np.array([[cos_transform, 0, sin_transform], 
                              [0, 1, 0], 
                              [-sin_transform, 0, cos_transform]])
    
    #We draw a cross at the origin:
    for i in range(-15, 15, 2):
        prod_3 = transform_arr @ np.array([0, 0, i])[:, np.newaxis]
        draw_point(prod_3.reshape(3) + centerOfMass, 3, outFile)
    
    #Then we draw the 2 planes, we are orthogonal to the best direction:
    for i in range(-20, 20, 5):
        point_3 = transform_arr @ np.array([i, 0, 0])[:, np.newaxis].reshape(3)
        draw_point(point_3 + centerOfMass, 3, outFile)
        point_4 = transform_arr @ np.array([0, i, 0])[:, np.newaxis].reshape(3)
        draw_point(point_4 + centerOfMass, 4, outFile)
    
        for j in range(-20, 20, 5):
            prod_1 = transform_arr @ np.array([i, j, dist_best_plane])[:, np.newaxis]
            prod_2 = transform_arr @ np.array([j, i, dist_best_plane - 15])[:, np.newaxis]
        
            point_1 = prod_1.reshape(3)
            point_2 = prod_2.reshape(3)
        
            draw_point(point_1 + centerOfMass, 1, outFile)
            draw_point(point_2 + centerOfMass, 2, outFile)
            



def calc_coord_plane(inputFile, best_direction, dist_best_plane, 
                     centerOfMass, precision):
    """To be able to have a nice PyMol output of our plan, we need to determine
    the coordinates of the 4 corners of our plane 
    (represented as a 3D rectangle)"""
    
    pdb_id = os.path.basename(inputFile).split('.')[0]
    outFile = open("./results/" + pdb_id + "_out.pdb", 'a')
    
    
    #Ici, on va utiliser la formule suivante:
    # cos(phi + pi/4) = R / H
    # <==> H = R / cos(phi + pi/4)
    
    #x = r * sin(phi) * cos(theta)
    #y = r * sin(phi) * sin(theta)
    #z = r * cos(phi)
    
    #Draw the axis orthogonal to the membrane
    #draw_axis(0, dist_best_plane, best_direction, centerOfMass, outFile)
    draw_axis(-30, 30, best_direction, centerOfMass, outFile)
        
    #draw_plane(dist_best_plane, idx_best_angles, precision, 
    #           centerOfMass, outFile)
    draw_plane(best_direction, dist_best_plane, centerOfMass, outFile)
    
    
    #We write a pymol_file, which can be executed automaticly:
    pymol_file = open("./results/cmd_pymol.pml", 'w')
    
    loading = "cmd.load('" + str(outFile.name) + "')\n"
    bg = "cmd.bg_color('white')\n"
    sticks = "cmd.show('sticks')\n"
    spheres = "cmd.show('spheres', 'HETATM')\n"
    
    for ligne in (loading, bg, sticks, spheres):
        pymol_file.write(ligne)
    
    outFile.close() ; pymol_file.close()
    
    
    
    
    
    
    #HETATM 2548  N   DUM  2548     -16.000  -6.000 -15.700
    #"{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(...)
    
    return

 



#MAIN
#if __name__ == "__main__":

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


# 4) Loop 
arr_nb_slices = np.zeros( (size_arr, size_arr) , dtype = int)




#Parallelization:
input = ( (i, j) for i, j in itr.product(range(size_arr), range(size_arr)) )
p = mp.Pool(mp.cpu_count())
results = p.map(parallelized_fun, input)
start_dist_plane_list = [result[0] for result in results]
mean_freq_list = [result[1] for result in results]
flat_start_dist_plane = np.array(start_dist_plane_list, dtype = int)
start_dist_plane_arr = flat_start_dist_plane.reshape( (size_arr, size_arr) )
arr_mean_freq = np.array(mean_freq_list).reshape( (size_arr, size_arr) )
p.close()
p.join()

#arr_mean_freq = np.zeros( (size_arr, size_arr) )
#for i in range(size_arr):
#    for j in range(size_arr):
        #if j > i+1:
#        arr_mean_freq[i, j] = parallelized_fun((i, j))




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


# 6) Finding the best position of the mb along the best direction:
best_start_dist_plane = start_dist_plane_arr[idx_theta, idx_phi]
print("START_DIST =", best_start_dist_plane)
idx_max_freq = improve_mb_position( best_direction, dict_CA, 
                                    nb_accessible_CA, list_resid )
dist_best_plane = best_start_dist_plane - idx_max_freq
print("DIST =", dist_best_plane)


# 6) Manage the output:

calc_coord_plane(inputFile, best_direction, dist_best_plane, 
                     centerOfMass, precision)

print("TOT_RUNTIME = ", time.time() - start_time, '\n')












