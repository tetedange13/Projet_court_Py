#!/usr/bin/env python3

#import math as m
import numpy as np

def polar_to_cartes(r, theta, phi):
    """Permet de passer des coord polaires aux coord cartesiennes (x;y;z), a 
    partir du rayon et des angles phi et theta"""
    return np.array([ r*np.sin(phi)*np.cos(theta), 
             r*np.sin(phi)*np.sin(theta),
             r*np.cos(phi) ])

def disque(pas):
    """Pour creer le disque, theta commence par valoir 0 
    et on fait varier phi. Prend en argument un pas, qui correspond un peu à une
    resolution (plus c'est grand, plus y a de vecteur"""
    
    for phi in np.arange(0.0, 2*np.pi-np.pi/pas-1, np.pi/pas):
        print(phi)
        
        
#MAIN

pas = 5
r = 3
phi = np.pi/2
resultats = np.zeros((pas, 2))


i=0
for theta in np.arange(0.0, np.pi/2, np.pi/pas):
    resultats[i, 0] = np.cos(theta)
    resultats[i, 1] = np.sin(theta)
    i += 1

    
print(resultats) #Afficher la facade
