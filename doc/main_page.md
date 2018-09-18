/*! \mainpage Welcome README Page

\section intro_sec Introduction
[doi_link]:https://doi.org/10.1093/bioinformatics/bth340
Welcome to the README Page of my Short Python Project.
The project consists in reimplementing an existing algorithm that finds
transmembrane parts from transmembrane proteins.
All the theorical concepts and the description of the method where found in the
following article:
> Gábor E. Tusnády, Zsuzsanna Dosztányi, István Simon; Transmembrane proteins 
> in the Protein Data Bank: identification and classification, Bioinformatics,
> Volume 20, Issue 17, 22 November 2004, Pages 2964–2972, [doi link][doi_link]
\n


\section install_sec Installation
\subsection step1 Step 1: Install all required Python packages
This program and all the associated modules are written in **Python 3**.\n
To work, this program needs the lastest following Python packages:
* multiprocess \version 0.70.6.1
* numpy \version 1.15.1
* biopython \version 1.72
> *REMARK:* The versions of the packages given are just INDICATIVE. The program 
> can as well work with later or previous versions \n
> (but some updates can sometimes be needed). \n
\n

\subsection step2 Step 2: Install Naccess
[naccess_link]: http://wolf.bms.umist.ac.uk/naccess/
The program needs Naccess software, in order to calculate the accessible surface
area (ASA). Sources for this software can be found on the [Naccess website]
[naccess_link].\n
(the archives are crypted, the password is "nac97", but don't tell anyone...) \n
\n

\subsection step3 Step 3: Get sources
The sources  can be found on my Git repository:
\verbatim 
git clone https://github.com/tetedange13/Projet_court_Py 
\endverbatim
But it's a bit messy... The best is still to use this nice archive that I 
dropped off on Moodle!\n
\n


\section use Usage
\subsection template Template of execution command line
The programm is called by executing "main.py". It takes several arguments, you
can have more information about the general use by typing
\verbatim
./main.py --[h]elp
\endverbatim \n

\subsection eg An example
Inside of the data/ folder, you have a couple of pdb files, on which you can try
the program. 1uaz and 6b87 are both TM protein, contrary to 1uw3, which is a
globular one. \n
Let's see an example:
\verbatim
./main.py -i ./data/6b87.pdb -p --naccess ../../Nacces/naccess --ASA 25
\endverbatim
Here we have set the precision ("-p") at 15, which will generate about 15*15
different directions of research. It should take around 3 secondes on a 2-cores
machine.
Depending on your own machine and the time you have, you can set a higher 
precision. \n
The script generate a **out pdb file, containing supplementary "DUM" atoms**, to
 represents the results of its search\n
It also generates a **PyMOL .pml file**, which contains a few commands for a
 nicelydisplay of the output pdb file. It can be launched like this:
\verbatim
pymol ./results/cmd_pymol.pml
\endverbatim


\n

\section troubles Troubleshooting
\subsection known Known Problems
Sometimes, it seems that the biopython module "Bio.NACCESS" have difficulties to
run Nacces by simply executing the "naccess" command...
So probably you will have to specify the path (absolute or relative) towards
the Naccess executable file.
Don't worry, it's very simple thanks to the --naccess option!

\subsection bug_report How to report a bug
If you have discovered a new bug (a one that I haven't already highlighted in my
written report or oral presentation), you can send me an email
<mailto:felix.deslacs@gmail.com>
  

