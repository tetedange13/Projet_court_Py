#!/usr/bin/env python3


"""
This module contains all the functions that deals with "drawing", i.e. appending
"DUM" HETEROATOMS in a pdb file (in the standard format of odb files), in order
to be able to have a nice output and maybe vizualise the execution 
"""



import os
import numpy as np



def draw_point(point, num, outFile):
    """
    Draw a point given its 3D coordinates and an int to identify the line
    
        Args: * A point <np.array> of size(1, 3)
              * A num <int> written in the file, to make the difference between 
                the atoms
              * An outFile <fileObject>, opened in appending, to write inside  
    
        Returns: None, it just appends a line to the outFile, corresponding to
                 the new HETATM DUM
    """
    
    well_formated = "{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}\n"

    outFile.write( well_formated.format( "HETATM", num, "N", "DUM", num,
                                         point[0],
                                         point[1],
                                         point[2] ))
    
    

def draw_axis(start, end, best_direction, outFile):
    """
    Draws an axis, by drawing a certain number of points (from start to end,
    by step of 5), along a direction given by a unit vector, by writting the
    coordinates in the outFile
    
    
    Args: * start <int> 
          * end <int>
          * best_direction <np.array> of size(1, 3)
     
    Returns: None, it just appends lines to the outFile, each corresponding to
             a new HETATM DUM
    """ 
    
    step = 5
    if start == 0 and end < 0: step *= -1    
                  
    for r in range(start, end, step):
        point_to_write = r * best_direction
        draw_point(point_to_write, r, outFile)

        

def draw_planes(best_direction, dist_best_plane, outFile):
    """Calculate the coordinates of the points belonging to the planes which
    are normal to the best direction and draw them inside the outFile
    
    
    Args: * best_direction <np.array> of size(1, 3)
          * dist_best_plane <float> distance between the best plane
            and the origin
          * An outFile <fileObject>, opened in appending, to write inside
    
    Returns: Returns: None, it is just writing inside a file
    
    
    REMARK: To manage this hard geometric task, we used a transformation matrix
            (transform_arr). It allow us to express simply the coordinates of
            a plane into the original (x,y,z) coordinates system of the pdb
            file. Then the application of the transformation matrix 
            transposes these coordinates into a new coordinate system, where
            the z axis corresponds to the axis normal to the membrane.
            
            The cos-sin of the angle between the best direction and the z axis
            have been calculated using both definition of the scalar production,
            as well as one cos-sin property (so a bit optimized): 
                    
                    cos_square(angle) + sin_square(angle) = 1
    """
    
    best_vect = np.array([best_direction]) * dist_best_plane
    produit = np.array([ [0.0, 0.0, 1.0] ]) @ best_vect.T
    cos_transform = produit[0, 0] / dist_best_plane                  
    sin_transform = np.sqrt(1 - cos_transform**2)
    transform_arr = np.array([[cos_transform, 0, sin_transform], 
                              [0, 1, 0], 
                              [-sin_transform, 0, cos_transform]])
    
    #We draw a cross at the origin:
    for i in range(-15, 15, 2):
        prod_3 = transform_arr @ np.array([ [0, 0, i] ]).T
        draw_point(prod_3.reshape(3), 3, outFile)
    
    #Then we draw the 2 planes, we are orthogonal to the best direction:
    for i in range(-20, 20, 5):
        point_3 = transform_arr @ np.array([ [i, 0, 0] ]).T.reshape(3)
        draw_point(point_3, 3, outFile)
        point_4 = transform_arr @ np.array([ [0, i, 0] ]).T.reshape(3)
        draw_point(point_4, 4, outFile)
    
        for j in range(-20, 20, 5):
            prod_1 = transform_arr @ np.array([ [i, j, dist_best_plane] ]).T
            prod_2 = transform_arr @ np.array([ [j, i, dist_best_plane - 15] ]).T
        
            point_1 = prod_1.reshape(3)
            point_2 = prod_2.reshape(3)
        
            draw_point(point_1, 1, outFile)
            draw_point(point_2, 2, outFile)

    
    print( "ANGLE = ", np.arccos(cos_transform) )
    
    

def esthetic_output(pdb_id, best_direction, dist_best_plane, precision):
    """
    Functions that executes all the "drawing" functions, in order to have a
    beautiful pdb_out
    
    It also generates a PyMOL file .pml, containing a success of commands to
    nicely print the pdb_out inside PyMOL
    
    
    Args: * pdb_id <str> (at least 4 characters)
          * best_direction <np.array> of size(1, 3)
          * dist_best_plane <float> distance between the best plane
            and the origin
          * precision <int>, approximately the number of directions explored,
            so can be considered as a resolution  
    
    Returns: None, as said before    
    """
    
    outFile = open("./results/" + pdb_id + "_out.pdb", 'a')
    
    #We draw the axis orthogonal to the membrane:
    draw_axis(-30, 30, best_direction, outFile)
    
    #Then we draw the 2 planes flanking the membrane:    
    draw_planes(best_direction, dist_best_plane, outFile)
    
    
    #We write a pymol_file, which can be executed automaticly:
    #pymol_file = open("./results/cmd_pymol_" + pdb_id + ".pml", 'w')
    pymol_file = open("./results/cmd_pymol.pml", 'w')
    
    loading = "cmd.load('" + outFile.name + "')\n"
    bg = "cmd.bg_color('white')\n"
    sticks = "cmd.show('sticks')\n"
    spheres = "cmd.show('spheres', 'HETATM')\n"
    
    for ligne in (loading, bg, sticks, spheres):
        pymol_file.write(ligne)
    
    outFile.close() ; pymol_file.close()
    
    
                
