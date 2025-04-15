mol load psf sin+dna_ions.psf pdb sin+dna_ions.pdb
set all [atomselect top all]
$all set beta 0.0
set sel [atomselect top "resname SIN"]
$sel set beta 1.0
set surf [atomselect top "resname SIN and \
((name \"SI.*\" and numbonds<=3) or (name \"N.*\" and numbonds<=2))"]
$surf set beta 10.0
$all writepdb sin+dna_restrain.pdb