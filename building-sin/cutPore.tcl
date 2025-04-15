# Author: Jeff Comer <jcomer2@illinois.edu>
# Cut a double-cone pore.

mol delete all
# Parameters:
set psf sin_bonded.psf
set pdb sin_bonded.pdb
set outPrefix sin_pore
set poreDiameter 24.0
set poreAngle 10.0

set pi [expr {4.0*atan(1.0)}]
set s0 [expr {0.5*$poreDiameter}]
set slope [expr {tan($poreAngle*$pi/180.0)}]

mol load psf $psf pdb $pdb
set pore [atomselect top "sqrt(x^2 + y^2) < $s0 + $slope*abs(z)"]
set poreAtoms [$pore get {segname resid name}]
$pore delete
mol delete top

# Use psfgen to delete the atoms.
package require psfgen
resetpsf
readpsf $psf
coordpdb $pdb
foreach atom $poreAtoms { delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2] }
writepsf $outPrefix.psf
writepdb $outPrefix.pdb

# Load the result.
mol load psf $outPrefix.psf pdb $outPrefix.pdb
