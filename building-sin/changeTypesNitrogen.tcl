set nit [atomselect top "type N"]
$nit set type NSI
set all [atomselect top all]
$all writepsf sin_pore_types.psf
$all writepdb sin_pore_types.pdb
