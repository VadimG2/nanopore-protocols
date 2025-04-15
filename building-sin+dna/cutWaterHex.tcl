# Author: Jeff Comer <jcomer2@illinois.edu>
# Cut the water to a hexagon.

# Parameters:
set sys sin+dna
set waterText water
set membraneText "resname SIN"

mol load psf ${sys}_sol.psf pdb ${sys}_sol.pdb
set sin [atomselect top $membraneText]
set minmax [measure minmax $sin]
set r [expr {0.5*([lindex $minmax 1 1]-[lindex $minmax 0 1])+2.0}]
set sqrt3 [expr {sqrt(3.0)}]
set cutText "($waterText) and ((abs(y) < 0.5*$r and abs(x) > 0.5*$sqrt3*$r) or (x > $sqrt3*(y+$r) or x < $sqrt3*(y-$r) or x > $sqrt3*($r-y) or x < $sqrt3*(-y-$r)))"
set cutSel [atomselect top $cutText]
set cutAtoms [lsort -unique [$cutSel get {segname resid}]]

package require psfgen
resetpsf
readpsf ${sys}_sol.psf
coordpdb ${sys}_sol.pdb
foreach atom $cutAtoms { delatom [lindex $atom 0] [lindex $atom 1] }
writepsf ${sys}_hex.psf
writepdb ${sys}_hex.pdb

mol delete all
mol load psf ${sys}_hex.psf pdb ${sys}_hex.pdb
