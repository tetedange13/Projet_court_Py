#!/usr/bin/env python3

from src.parse_files import *
from src.geometry import *
from src.drawing import *

import sys
import time
from getopt import getopt
import itertools as itr
import multiprocess as mp



#To simplify the syntax, in the next functions, the use of "CA" in any variable
#will refer to "all the CA which are enough accessible to the solvent":

    
def start_position_plane(unit_vect, dict_CA, nb_CA, list_resid):
    """
    Takes an array containing the coord of a unit vect.
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
    will be used to slide along the axis (function next next)
    
    Args: * A unit vector i.e.a direction in the 3D space (<np.array> type,
            with size(1, 3))
          * The nested dict with information about the CA (+ 3D coordinates) 
          * The number of ACCESSIBLE CA (we do not precise anymore 
            at this stage) 
          * The list of the resid (already explained, as <list> type)   
    
    Returns: * The number of slices i.e. the number of steps needed to cross
               the protein from a side to the other, along the best axis given
               (as a <float>)
             * The array with starting point distance array, which will be
               decremented to simulate the sliding
             * The dist between the plane in its starting position and the
               origin of the system (will be used later in the improve
               function)      
    """

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
    """
    Reads the dist_arr, in order to find which one are in the current slice.
    Then by the resid associated to the index (via list of resid), goes into
    the Naccess output file, to know if the corresponding amino acid is 
    hydrophobic or not
   
    
    Args: * The array with starting point distance array, which will be
            examinated to look for values which are in the intervall of
            [0, 10] A  (type <np.array>)
          * The number of ACCESSIBLE CA (we do not precise anymore 
            at this stage)
          * The nested dict with information about the CA (+ 3D coordinates)
          * The list of the resid (already explained, as <list> type)   
          
    Returns: Frequence of hydrophobic among all accessible aa 
             in the current slice
    """

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
    """
    Takes a slice positionned at its start point, by taking 
    its start_dist_points_arr. 
    Slides this slice nb_slices times, by decrementing the dist_arr of 1A.
    At each step, the frequency of hydrophobic aa among the accessible aa is
    calculated, in order to have the mean of the frequency of hydrophobic aa
    along the current direction
    
    Args: * The array with starting point distance array, which will be
            decremented to simulate the sliding
          * The number of ACCESSIBLE CA (we do not precise anymore 
            at this stage)
          * The nested dict with information about the CA (+ 3D coordinates) 
          * The list of the resid (already explained, as <list> type)  
    
    Returns: The mean of hydrophobic frequency, weighted by the square of the
             number of slices, in order to take in count the big differences of
             size among all the different number of slices, which used to bias
             the comparaison of the means
    
    REMARK: The different values of frequency along a given direct are not kept
    individually, in order to reduce the cost in memory.
    That is why, after found the good normal vect to the membrane, we will have 
    to slide one last time along this direction, in order to well position the
    membrane
    """

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
    """
    To calculate the mean freq_hydrophob on all the directions 
    (in a parallelized way)
    
    Args: A tuple of indexes which will give a direction, on which a mean of
          hydrophobic frequency will be calculated
    
    Returns: * The initial distance between the plane (positionned near from the
               closest point along this direction) and the origin of the system
             * The value of the mean of hydrophobic frequency along 
               this direction 
    """
    
    i,j = tupl
    unit_vect = arr_unit_vect[:, i, j]
    
            
    #4-1) We position the plan well:
    tuple_returned = start_position_plane( unit_vect, dict_CA, nb_accessible_CA,
                                           list_resid )
    nb_slices, start_dist_points_arr, start_dist_plane = tuple_returned
        
    #4-2) We slide along the current direction:   
    return ( start_dist_plane, 
             sliding_slice(nb_slices, start_dist_points_arr, 
                           nb_accessible_CA, dict_CA, list_resid) )         



def improve_mb_position(best_direction, dict_CA, nb_CA, list_resid):
    #Je l'ai un peu optimisée, pour éviter de stocker en memoire une matrice de
    #freq, que je dois ensuite parcourir pour trouver le max avec np.argmax()
    #Par contre, le repositionnement du plan au debut, il est pas tres 
    #optimal, puisque que tu dois recalculer la dist entre le plan et tous les 
    #points...

    """
    As we did not keep in memory all the values of frequencies along all the
    different direct, we have to slide along the best direction one more time,
    in order to find the optimal position of the membrane along the best vector
    
    Aim si to find the position where the associated slice has the maximum
    hydrophobic frequency
    
    
    Args: * best_direction <np.array> of size(1, 3)
          * The nested dict with information about the CA (+ 3D coordinates)
          * The number of ACCESSIBLE CA (we do not precise anymore at this stage)
          * The list of the resid (already explained, as <list> type)  
    
    Returns: The index of the position of the membrane along the best axis, 
             which will be used to calculate the actual distance between 
             the membrane and the origin
    """
    
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
    
    


#MAIN
if __name__ == "__main__":

    start_time = time.time()

    # 0) We parse the arguments given:
    path_to_naccess_exe = 'naccess'
    str_usage = "Usage: main.py -i <inputFile.pdb> --[n]access" \
                " <path_to_naccess_exe> --[p]recision <int>" \
                " --[a]sa <float(thresoldASA)>"

    opts, args = getopt(sys.argv[1:], "hani:p:", [ "naccess=", 
                                                   "precision=",
                                                   "asa=",
                                                   "help="])

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(str_usage)
            sys.exit()
            
        elif opt == "-i":
            inputFile = arg  
            
        elif opt in ("-n", "--naccess"):
            path_to_naccess_exe = arg
        
        elif opt in ("-p", "--precision"):
            precision = int(arg)
            
        elif opt in ("-a", "--asa"):
            threshold_ASA = float(arg)    
            
        else:
            sys.exit("ERROR: Invalid argument passed\n", str_usage)        


    # 1) We extract the needed data:
    dict_CA, nb_accessible_CA, nb_tot_CA = get_acessible_CA(inputFile, 
                                                            path_to_naccess_exe)

    dict_CA, list_resid, centerOfMass = get_coord( dict_CA, nb_accessible_CA, 
                                                   nb_tot_CA, inputFile )


    # 2) We transform the coordinates to set the center of mass as the origin;
    #As well as in the dict as in the output file:

    generate_transformed_outFile(inputFile, centerOfMass)
    dict_CA = transform_coord(dict_CA, centerOfMass)


    # 3) We generate the  arrays (only once), using these formulas
    #To swith from polar to cartesian coordinates:
    #x = r * sin(phi) * cos(theta)
    #y = r * sin(phi) * sin(theta)
    #z = r * cos(phi)

    arr_cos_theta, arr_sin_theta = generate_trigo_arr(precision)
    arr_cos_phi, arr_sin_phi = arr_cos_theta, arr_sin_theta

    #Use of matricial product to calculate each of the x, y and z matrixes
    #of the unit vector dividing half of the the 3D space
    #(necessity to reshape the different matrixes, to make the product)
    size_arr = precision + 1
    vectArr_X = arr_cos_theta.T @ arr_sin_phi
    vectArr_Y = arr_sin_theta.T @ arr_sin_phi
    vectArr_Z = np.tile( arr_cos_phi, (size_arr, 1) )

    #REMARK: The vectArr_Z has been created by repeating the line 
    #of cos(phi) as many times there are values for the theta angle

    #We merge the 3 arrays in a single one, to simplify its passing 
    #as an argument: 
    arr_unit_vect = np.zeros( (3, size_arr, size_arr), dtype = float )
    arr_unit_vect[0, :, :] = vectArr_X
    arr_unit_vect[1, :, :] = vectArr_Y
    arr_unit_vect[2, :, :] = vectArr_Z


    # 4) Parallelized nest for loop, which :
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
    idx_max_freq = improve_mb_position( best_direction, dict_CA, 
                                        nb_accessible_CA, list_resid )
    dist_best_plane = best_start_dist_plane - idx_max_freq


    # 7) Manage the output:
    pdb_id = os.path.basename(inputFile).split('.')[0]
    esthetic_output(pdb_id, best_direction, dist_best_plane, precision)

    print("TOT_RUNTIME = ", time.time() - start_time, '\n')












