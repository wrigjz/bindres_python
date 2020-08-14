reint
# Load the pdb files and get the binding residues
load post_minix.pdb
# Load speedfill pdb and match critires w speedfill w binding residues
load speedfill.pdb
cmd.remove("(all) and hydro")
select SPH, /speedfill///SPH/
select wanted, byres post_minix within 5.0 of SPH
select wanted_ca, wanted and name ca
iterate wanted_ca, print(resi, resn, chain)
