#!!! {0 4 65} FOR NEGATIVE, {0 4 -115} for POSTIVE/NEUTRAL!!! 
mol load psf sin+dna.psf pdb sin+dna.pdb
set all [atomselect top all]
set sel [atomselect top "segname P1"]
$sel moveby {0 4 65}
$all writepdb sin+dna_placed.pdb
