#!/usr/bin/env python3

#import math as m
import numpy as np
import sys

def polar_to_cartes(r, theta, phi):
    """Permet de stepser des coord polaires aux coord cartesiennes (x;y;z), a 
    partir du rayon et des angles phi et theta"""
    return np.array([ r*np.sin(phi)*np.cos(theta), 
             r*np.sin(phi)*np.sin(theta),
             r*np.cos(phi) ])
             

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
        
        
        
#MAIN

#Formulas, from polar to cartesian:
#x = r * sin(phi) * cos(theta)
#y = r * sin(phi) * sin(theta)
#z = r * cos(phi)

r = 3
step = int(sys.argv[1])



arr_cos_theta, arr_sin_theta = rempli_arr(step, "theta")
arr_cos_phi, arr_sin_phi = rempli_arr(step, "phi")


#Use of matricial product to calculate each of the x, y and z matrixes
#(necessity to reshape the phi matrixes, to make the product)
size_arr_phi = np.shape(arr_cos_phi)[0]
size_arr_theta = np.shape(arr_cos_theta)[0]  
X_arr = r * arr_cos_theta @ arr_sin_phi.reshape(1, size_arr_phi)
Y_arr = r * arr_sin_theta @ arr_sin_phi.reshape(1, size_arr_phi)
Z_arr = r * np.array( [arr_cos_phi.reshape(1, size_arr_phi)] * size_arr_theta )

#REMARK: The Z_arr has been created by repeating the line of cos(phi) as many
#times there are values for the theta angle
    
print(X_arr)
print("\n ESPACE\n")
print(Y_arr)
print("\n ESPACE\n")
print(Z_arr)












