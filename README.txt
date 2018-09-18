===========================
   Short Python Projects
   
Topic: Detection of 
transmembrane regions

AUTHOR: Felix Vandermeeren
YEAR: (sept) 2018-2019
===========================


================
EMERGENCY README
================


Emergency description
---------------------
The project consists in reimplementing an existing algorithm that finds
transmembrane parts from transmembrane proteins
From this article: https://doi.org/10.1093/bioinformatics/bth340


Emergency procedure
-------------------

Normally, you should have the fancy documentation, generated with Doxygen,
by typing in a terminal from the root directory of the project:

            firefox ./doc/html/index.html
            
            
If this command does not work, and that you have Doxygen installed on your
machine, you can regenerate the documentation with:

            doxygen doxy.conf  
And then:   firefox ./doc/html/index.html 


If it is still not working, the essential of the documentation is below.
Just imagine that it is nice and beautiful:


Requirements
------------
* Naccess program is needed (http://wolf.bms.umist.ac.uk/naccess/, psswd = "nac97")

* As for the following python modules: numpy, biopython and multiprocess
(not a precise version needed)



Emergency content
-----------------

Help can be found with:
    ./main.py --[h]elp
    
    >> "Usage: main.py -i <inputFile.pdb> --[n]access" \
                " <path_to_naccess_exe> --[p]recision <int>" \
                " --[a]sa <float(thresoldASA)>"
                

Example of how to run the script, with all the parameters:
    ./main.py -i ./data/6b87.pdb -p 50 --naccess ../../Nacces/naccess --ASA 25



