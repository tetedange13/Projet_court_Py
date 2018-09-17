#!/usr/bin/env python3

import numpy as np


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


