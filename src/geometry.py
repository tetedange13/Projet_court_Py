#!/usr/bin/env python3


"""
This module contains all functions related to the manipulations of geometrical 
objects like points (can be numpy arrays or dict), vectors (more often numpy
arrays), or planes.
"""



import numpy as np



def generate_trigo_arr(precision):
    """
    Share the half of the trigonometric circle in (precision+1) angles and
    generates 2 arrays, corresponding to the sin and cos of these angles
     
     
    Args: precision <int>, approximately the number of directions explored,
          so can be considered as a resolution  
     
    Returns: A tuple of 2 <np.arrays>, both of size(1, (precision+1))
    
    
    REMARK: The function is a bit optimised, on the fact hat we calculate a cos
    and sin only for 1 dial on 2.
    The other being deduced from the cos-sin properties as follows:
    
        cos(pi-x) = -cos(x) ; sin(pi-x) = sin(x) <=> N-O dial
    """
    
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
        # "-1/precision" found a bit empirically..
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
    
    return ( np.array([arr_cos]), np.array([arr_sin]) )        



def dist_point_plane(point, plane):
    """Calculates the distance between a plane of form (a, b; c; d) and
    a point of form (x; y; z)
    
    
    Args: * point <dict> wih keys = ['x', 'y', 'z']
          * plance <np.array> of size(1, 4)
          
    Returns: The distance as a float
    
    
    REMARK: As we are working with planes moving along unit vector with norm R,
            the d coefficient of any plane -R_square
    """
    
    arr_point = np.hstack( (point['x'], point['y'], point['z'], 1) )
        
    numer = abs( arr_point @ plane.T )
    denom = np.sqrt( plane[0:3, ] @ plane[0:3, ].T )
    return numer/denom


