#!/usr/bin/env python3

#import math as m
import numpy as np
import sys

def get_C_alpha_coord(pdbFile):
    """Get the coordinate of the C-alpha from a pdb file (file obj) 
    given as argument. Count also the number of C-alpha"""
    
    coord_dict = {}
    nb_Calpha = 0
    for line in pdbFile:
        if line.startswith( "ATOM" ) and line[12:16].strip() == "CA":
            nb_Calpha += 1
            resid = int( line[22:26].strip() )
            coord = { 'x': float( line[30:38].strip() ), 
                      'y': float( line[38:46].strip() ), 
                      'z': float( line[46:54].strip() ) 
                    }
            coord_dict[resid] = coord
    
    return (coord_dict, nb_Calpha-1)


def calc_centerOfmass(dict_coord, nb_Calpha):
    """Calculate the centroid (center of mass) of a protein, given the 
    coordinates of all its C-alpha, and return a transformed coord dict"""
    
    x_sum, y_sum, z_sum = 0, 0, 0 
    for resid in dict_coord.keys():
        x_sum += dict_coord[resid]['x']
        y_sum += dict_coord[resid]['y']
        z_sum += dict_coord[resid]['z']
    
    centerOfMass = np.array( [x_sum, y_sum, z_sum]) / nb_Calpha
    
    tranformed_dict = dict_coord
    for resid in dict_coord.keys():
        tranformed_dict[resid]['x'] -= centerOfMass[0]
        tranformed_dict[resid]['y'] -= centerOfMass[1]
        tranformed_dict[resid]['z'] -= centerOfMass[2]
    
    return (centerOfMass, tranformed_dict)
    
             

def rempli_arr(step, angle):
    """Prend un pas (sorte de resolution), decoupe le cercle trigo et renvoie 2
     matrices colonnes des cos et sin de ces angles
     Fonctionne aussi bien pour theta que pour phi"""
     
    #Prop of trigo functions trigo used:
    #cos(pi-x) = -cos(x) ; sin(pi-x) = sin(x) <=> N-O dial
    #cos(pi+x) = -cos(x) ; sin(pi+x) = -sin(x) <=> S-O dial
    #cos(-x) = cos(x) ; sin(-x) = -sin(x) <=> S-E dial
    
    
    #We 1st disinguish the cases where we have to deal with theta (from 0 to 2pi)
    #Or phi (from à to pi):
    if angle == "theta":
        arr_cos = np.zeros((2*step, 1), dtype=float)
        arr_sin = np.zeros((2*step, 1), dtype=float)
        nb_dials = 4
        
    elif angle == "phi":       
        arr_cos = np.zeros((step+1, 1), dtype=float)
        arr_sin = np.zeros((step+1, 1), dtype=float)
        nb_dials = 2 
    
    else:
        sys.exit("Pas bon argument donne pour l'angle")
    

    #We start with filling with evident values (not zero) of cos et sin:
    #En 0:
    arr_cos[0, 0] = 1.0
    #En pi:
    arr_cos[step, 0] = -1.0
    
    if step%2 == 0: #If the step is even, there are 2 more (pi/2 and 3pi/2)
        idx_pi_over_2 = int( (step+1)/2 )
        #En pi/2:
        arr_sin[idx_pi_over_2, 0] = 1.0
        if angle == "theta": #This value won'be needed for phi
            #En 3pi/2:
            arr_sin[3*idx_pi_over_2, 0] = -1.0        
            
            
    #N-E DIAL (will serve as a reference to fill the other dials):
    i=1
    for angle in np.arange(np.pi/step, np.pi/2 - 1/step, np.pi/step):
    #J'ai trouve le "-1/pas" un peu empiriquement => A VOIR...
        arr_cos[i, 0] = np.cos(angle)
        arr_sin[i, 0] = np.sin(angle)
        i += 1     
    
    nb_to_calc = i - 1 #Number of values to calculate
    
    
    #Then deduction of the other dials, based on the N-E one:
    for dial_idx in range(1, nb_dials): #Starts at 1, cuz NE dial has index 0
        sign_cos = 1 ; sign_sin = 1 #Initialized at the beginning of each turn
    
        #Depending on the parity of the step, we don't start from the same index
        if step%2 == 0:  
            start = dial_idx * idx_pi_over_2 + 1
        else:
            start = dial_idx * i        
        
        if dial_idx == 2: #No need to reverse in the S-O dial
            sign_cos = -1; sign_sin = -1
            end = start + nb_to_calc
            arr_cos[start:end, 0] = sign_cos * arr_cos[1:i, 0]
            arr_sin[start:end, 0] = sign_sin * arr_sin[1:i, 0]        
        
        else: #Reversion for the 2 other dials, to reflect trigo circle
            if dial_idx == 1: #N-O DIAL
                sign_cos = -1

            elif dial_idx == 3: #S-E DIAL
                sign_sin = -1
                if step%2 != 0: #Particular case when S-E dial and odd step
                    start = dial_idx*i - 1  
        
            end = start + nb_to_calc
            arr_cos[start:end, 0] = sign_cos * arr_cos[1:i, 0][::-1]
            arr_sin[start:end, 0] = sign_sin * arr_sin[1:i, 0][::-1]

    
    return (arr_cos, arr_sin)        


#def plane_equation():
    
        

def dist_point_plane(point, plane):
    """Prend un point provenant du pdb (donc un dict)
    et un plan = (a, b; c; d), avec:
        -- (a; b; c)  = (x; y; z) du point au bout d'un vecteur 
        -- et d = rayon de la sphère"""
    
    arr_point = (x0; y0; z0; 1)
        
    numer = abs( point @ plane.T )
    denom = np.sqrt( plane[0:4, ] @ plane[0:4, ].T )
    return numer/denom
    

def is_in_slice(point, planeDown, planeUp):
    """Determine if a given C-alpha (aka point with (x;y;z) coord) is inside a
    given slice of 1A wide (between planeI and planeI_plus1)"""
    
    if dist_point_plane(point, planeDown) < 1 \
                                and dist_point_plane(point, planeUp) < 1:
        return 1    
    
    else:
        return 0                                
    
def fonction (r):
    size_arr_theta = np.shape(arr_cos_theta)[0]  
    X_arr_down, X_arr_up = r * arr_cos_theta @ arr_sin_phi.T, \
                           (r+1) * arr_cos_theta @ arr_sin_phi.T
    Y_arr_down, Y_arr_up = r * arr_sin_theta @ arr_sin_phi.T, \
                           (r+1) * arr_sin_theta @ arr_sin_phi.T
    Z_arr_down, Z_arr_up = r * np.tile( arr_cos_phi.T, (size_arr_theta, 1) ), \
                           (r+1) * np.tile( arr_cos_phi.T, (size_arr_theta, 1) )  
    
    point = np.array([)
    taille_i, taille_j = np.shape(X_arr_down)
    for i in range(taille_i):
        for j in range(taille_j):
            planeDown = np.array( [X_arr_down[i,j],
                                   Y_arr_down[i,j],
                                   Z_arr_down[i,j],
                                   -r**2 ] )
            planeUp = np.array( [X_arr_up[i,j],
                                   Y_arr_up[i,j],
                                   Z_arr_up[i,j],
                                   -r**2 ] )
            print(is_in_slice(
         
#MAIN

#Formulas, from polar to cartesian:
#x = r * sin(phi) * cos(theta)
#y = r * sin(phi) * sin(theta)
#z = r * cos(phi)

r = 3
step = int(sys.argv[1])


#These arrays need to be created only one:
arr_cos_theta, arr_sin_theta = rempli_arr(step, "theta")
arr_cos_phi, arr_sin_phi = rempli_arr(step, "phi")


#Use of matricial product to calculate each of the x, y and z matrixes
#(necessity to reshape the phi matrixes, to make the product)
size_arr_theta = np.shape(arr_cos_theta)[0]  
X_arr = r * arr_cos_theta @ arr_sin_phi.T
Y_arr = r * arr_sin_theta @ arr_sin_phi.T
Z_arr = r * np.tile( arr_cos_phi.T, (size_arr_theta, 1) ) 

#REMARK: The Z_arr has been created by repeating the line of cos(phi) as many
#times there are values for the theta angle

mon_plan = np.array( [X_arr[2,2], Y_arr[2,1], Z_arr[2,1], -r**2] )
mon_point = np.array( [1, 2, 3, 1] ) 
#print( dist_point_plane(mon_point, mon_plan) )

mon_pdb = open("/home/sdv/m2bi/fvandermeeren/Bureau/6gx6.pdb", 'r')
coord_mon_pdb, nb_Calpha = get_C_alpha_coord(mon_pdb)
mon_pdb.close()


print(nb_Calpha)
print(len(coord_mon_pdb))

#for resid in coord_mon_pdb:
#    print(coord_mon_pdb[resid])

centerOfMass, transformed_coord = calc_centerOfmass(coord_mon_pdb, nb_Calpha)

print( calc_centerOfmass(transformed_coord, nb_Calpha)[0] )


fonction(3)
    
#print(X_arr)
#print("\n ESPACE\n")
#print(Y_arr)
#print("\n ESPACE\n")
#print(Z_arr)












