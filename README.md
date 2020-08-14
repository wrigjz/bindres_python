# bindres_python
###################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
###################################################################################################
This is a implimentation of a simple method to try to find which residues on a single PDB chain 
Are invovled in binding to abother protein - it requires only the free single chain of one of the 
chains in the complex

Requirements are AmberTools, hbplus, Speedfill and FreeSASA

the directory structure should be something like this

top/bindres_scripts

   /bindresscripts (If you want to use our consurf implimentation)

   /pdbidchain    e.g. 1fnfa

              /input.pdb   this is a single chain from a pdb file

the pdb file can be prepared using the get_pdb_get_from_archive.py script by giving it a file with a list of pdb and chain codes in a format such as 1TSR B etc etc, you may need to edit the script to point to your own local pdb archive


if you cd into pdbidchain and run
        ../bindres_scripts/bindres.sh

then bindres should run, the default is to also run a local version of consurf,  which can be ignored by editing the bindres.sh file to comment out a few lines and uncomment some others, you can also choose to run consurf and use the SEQRES records instead
