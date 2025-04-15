package require psfgen
resetpsf
readpsf ../building-dna/polys_autopsf.psf
coordpdb ../building-dna/polys_autopsf.pdb
readpsf ../building-sin/sin_pore_charges.psf
coordpdb ../building-sin/sin_pore_charges.pdb
writepsf sin+dna.psf
writepdb sin+dna.pdb
