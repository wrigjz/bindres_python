###################################################################################################
## Jon Wright, IBMS, Academia Sinica, Taipei, 11529, Taiwan
## These files are licensed under the GLP ver 3, essentially you have the right
## to copy, modify and distribute this script but all modifications must be offered
## back to the original authors
###################################################################################################
import math

# This script will take a pdb file prepared by AMBER and a speedfill pdb SPH file
# and return the RESIDues that are within a cutoff DISTANCE from any of the SPH 'atoms'
# the post_mini_noh.pdb file should already have had the hydrogens removed

# Set the cutoff limit
CUTOFF = 5.0

# Get a line count for the arrays
NUM_LINES = 1
NUM_LINES += sum(1 for lines in open("post_mini_noh.pdb", "r"))
SPH_LINES = 1
SPH_LINES += sum(1 for lines in open("speedfill.pdb", "r"))

PDBFILE = open("post_mini_noh.pdb", "r")
SPHFILE = open("speedfill.pdb", "r")

ATOMID = [0 for i in range(0, NUM_LINES)]
ATOMNA = [0 for i in range(0, NUM_LINES)]
RESNA = [0 for i in range(0, NUM_LINES)]
CHAIN = [0 for i in range(0, NUM_LINES)]
RESID = [0 for i in range(0, NUM_LINES)]
X_CRD = [0 for i in range(0, NUM_LINES)]
Y_CRD = [0 for i in range(0, NUM_LINES)]
Z_CRD = [0 for i in range(0, NUM_LINES)]
SPHID = [0 for i in range(0, SPH_LINES)]
X_SPH = [0 for i in range(0, SPH_LINES)]
Y_SPH = [0 for i in range(0, SPH_LINES)]
Z_SPH = [0 for i in range(0, SPH_LINES)]

# Set up some counters
INDEX = -1
RESCOUNT = -1
PREVIOUS_RESID = 99999

# Reas in the PDB file and store the cords, RESIDue id etc
for LINE in PDBFILE:
    if LINE[0:4] == "ATOM" or LINE[0:6] == "HETATM":
        INDEX += 1
        # Save the info we need for later
        ATOMID[INDEX] = LINE[6:11]
        ATOMNA[INDEX] = LINE[12:16]
        RESNA[INDEX] = LINE[17:20]
        CHAIN[INDEX] = LINE[21:22]
        RESID[INDEX] = int(LINE[22:27].strip())
        X_CRD[INDEX] = float(LINE[30:38])
        Y_CRD[INDEX] = float(LINE[38:46])
        Z_CRD[INDEX] = float(LINE[46:54])
        if RESID[INDEX] != PREVIOUS_RESID:
            RESCOUNT += 1
            PREVIOUS_RESID = RESID[INDEX]

# Once we have that then we can set up the arrays to keep the results
# Need to use +2 as the RESIDue count starts at 1
WANTED_RESID = [0 for i in range(0, RESCOUNT+2)]
WANTED_RESNA = [0 for i in range(0, RESCOUNT+2)]

# Read in the Speedfill file and store the crds and sphere number
for LINE in SPHFILE:
    if LINE[0:4] == "ATOM" or LINE[0:6] == "HETATM":
        X_SPH = float(LINE[30:38])
        Y_SPH = float(LINE[38:46])
        Z_SPH = float(LINE[46:54])
        # Now check atoms to see if there are within cutoff DISTANCE of a sphere
        for i in range(0, INDEX+1):
            if RESNA[i] != "ACE" and RESNA[i] != "NME": # Ignore NME and ACE
                X_DIFF_SQ = (X_CRD[i] - X_SPH) ** 2
                y_DIFF_SQ = (Y_CRD[i] - Y_SPH) ** 2
                Z_DIFF_SQ = (Z_CRD[i] - Z_SPH) ** 2
                DISTANCE = math.sqrt(X_DIFF_SQ + y_DIFF_SQ + Z_DIFF_SQ)
                if DISTANCE <= CUTOFF: # Flag this residue as wanted
                    #print(DISTANCE, ATOMID[i], SPHID[j])
                    WANTED_RESID[RESID[i]] = RESID[i]
                    WANTED_RESNA[RESID[i]] = RESNA[i]

# Now check the residues and print out the ones we need
for i in range(0, RESCOUNT+1):
    if WANTED_RESID[i] != 0:
        print(WANTED_RESID[i], WANTED_RESNA[i], "X")
