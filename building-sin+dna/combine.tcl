package require psfgen
resetpsf
readpsf ../building-dna/dsdna.psf
coordpdb ../building-dna/dsdna.pdb
readpsf ../building-sin/sin_pore_charges.psf
coordpdb ../building-sin/sin_pore_charges.pdb
writepsf sin+dna.psf
writepdb sin+dna.pdb
