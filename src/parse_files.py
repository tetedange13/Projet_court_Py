#!/usr/bin/env python3


from Bio.PDB import PDBParser 
from Bio.PDB.NACCESS import run_naccess


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


