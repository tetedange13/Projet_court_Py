#!/usr/bin/env python3


"""
This module contains all the functions implied in the parsing and processing
of the files which are needed to generate an output:

* Collection of the relative ASA for each residu from the .rsa file, in order to
determine which are the amino acids which can be considered as accessible to the
solvant (so like to the membrane);

* Collection of the 3D coordinates of these accessible CA, for further 
manipulations.
"""



import os
from Bio.PDB import PDBParser 
from Bio.PDB.NACCESS import run_naccess
import numpy as np



def get_acessible_CA(inputFile, path_to_naccess_exe, thresold_ASA=30):
    """
    Runs Naccess on a pdb file, by using the Bio.Naccess (sub)module from the
    Biopython module. 
    Parses the generated output, in order to get the relative ("REL") ASA on all
    the input protein.
    
    Given a threesold on the value of ASA, descrimines the residus between the
    accessible to the solvant and the others
    Only the accessible one are kept and some of their information (like ASA,
    chain, etc) are kept in a nested dictionnary 
    
    
    Args: * An input pdb file as form of <str> (path to the file)
          * The path to the Naccess executable file (as a <str>)
          * The threesold of ASA beyond which a residu can be considered as
            accessible to the solvent (optionnal <float>, default value at 30)
            
    Returns: * A nested dict with information about the CA 
               (no coordinates for now)
             * The number of accessible CA as <int>, will be used in the next
               function
             * The total number of CA as an <int> (also for the next function)  
    """
    
    #If NACCESS available:
    pdb_fileName = os.path.basename(inputFile)
    pdb_id = pdb_fileName.split('.')[0]
    #p = PDBParser(QUIET = True)
    #structure = p.get_structure(pdb_id, inputFile)
    #model = structure[0]
    #output_naccess = run_naccess(model, inputFile, 
    #                             naccess = path_to_naccess_exe)[0]
    
    #If NACCESS NOT available:
    results_naccess = os.system ("~/Naccess/naccess " + inputFile)
    os.system("rm -f ./" + pdb_id + ".asa")
    os.system("rm -f ./" + pdb_id + ".log")
    os.system("mv ./" + pdb_id + ".rsa ./data/")
    naccessFile = open("./data/" + pdb_id + ".rsa", 'r')
    output_naccess = naccessFile.readlines()
    naccessFile.close()
     
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
    """
    Updates the dictionnary containing the information about the ACCESSIBLE
    (only) CA, by appending their (x,y,z) coordinates, got from the input 
    pdb file
    
    By the way calculates the center of mass of the protein and generate a list
    of the resid of the CA, in order to later access rather by their index, to
    correspond avec the index of matrixes
    
    
    Args: * A nested dict with information about the CA
            (no coordinates for now)
          * The nb_accessible_CA as an <int>, to initialise the resid list
          * The total number of CA, needed for the calculation of the 
            center of mass (as an <int>)
          * An input pdb file, from which the (x,y,z) coordinates of the 
            accessible CA (only) will be extracted
    
    Returns: * A nested dict with information about the CA
               (coordinates added now)
             * The list of the resid (already explained, as <list> type)  
             * The coordinates of the center of mass, as a <np.array> of 
               size(1, 3)
    """

    dict_CA_final = dict_CA #To avoid side effects
    
    pdbFile = open(inputFile, 'r')
    pdb_id = os.path.basename(inputFile).split('.')[0]
    
    outFile = open("./results/" + pdb_id + '_out.pdb', 'w')
    
    #To index the resid
    list_resid = [0] * nb_accessible_CA ; i = 0
    #To calculate the center of mass of the protein
    sum_X, sum_Y, sum_Z = 0, 0, 0
    
    resid = 0 #Problems with resid starting at 3, or chains etc in pdbFile
    for line in pdbFile:
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
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


    pdbFile.close()
    centerOfMass = np.array([sum_X/nb_tot_CA, sum_Y/nb_tot_CA, sum_Z/nb_tot_CA]) 
    return (dict_CA_final, list_resid, centerOfMass)



def transform_coord(dict_coord, centerOfMass):
    """
    Transforms in situ the coordinates of CA inside a dict, in order to set the
    center of mass  as the new origin of the coordinate system.
    
    
    Args: * dict_coord <dict> whose keys are resid, pointing towards another
            nested dict containing several info and in particular the
            (x,y,z) coordinates of the CA
          * centerOfMass <np.array> of size(1, 3)
          
    Returns: A copy of the dict passed as argument, but with coordinates
             transformed to set the center of mass as origin 
    """
    
    transformed_dict = dict_coord #To avoid side effects
    for resid in dict_coord.keys():
        transformed_dict[resid]['x'] -= centerOfMass[0]
        transformed_dict[resid]['y'] -= centerOfMass[1]
        transformed_dict[resid]['z'] -= centerOfMass[2]
        
        
    return transformed_dict



def generate_transformed_outFile(inputFile, centerOfMass):
    """
    To get rid of questions like "have I thought to transform all the coordinates
    with the center of mass every time I write a point?"; or "Do I need to use 
    the transformed pr untransformed coordinates to manage this scalar product?";
    
    Thanks to this function, the coordinates of ALL the output pdb file are
    transformed in order to set the center of mass of the protein as the origin
    of the coordinates system
    
    
    Args: * inputFile <str> path to the file to take the original 
            coordinates from
          * centerOfMass <np.array> of size(1, 3) 
    
    Return: None, except a brand new transformed pdb output file
    """
    
    pdb_id = os.path.basename(inputFile).split('.')[0]
    pdbFile = open(inputFile, 'r')
    outFile = open("./results/" + pdb_id + '_out.pdb', 'r+')
    tmp_dict = {}
    arbitrary_number = 12

    for line in pdbFile:
        if line.startswith("ATOM"):
            X_coord = float( line[30:38].strip() )
            Y_coord = float( line[38:46].strip() )
            Z_coord = float( line[46:54].strip() )
            
            tmp_dict[arbitrary_number] = { 'x': X_coord,
                                           'y': Y_coord,
                                           'z': Z_coord }
            tmp_dict = transform_coord(tmp_dict, centerOfMass)

            to_add = "{:8.3f}" \
                     "{:8.3f}{:8.3f}".format( tmp_dict[arbitrary_number]['x'],
                                              tmp_dict[arbitrary_number]['y'],
                                              tmp_dict[arbitrary_number]['z'] ) 
            newLine = line[:30] + to_add + line[54:]
            outFile.write(newLine)
            
    outFile.write("MASTER\n")
    pdbFile.close() ; outFile.close()



