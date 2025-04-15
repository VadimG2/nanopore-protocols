# Author: Jeff Comer <jcomer2@illinois.edu>
# Add ions (KCl) of given ionic strength.
# vmd -dispdev text -e addIons.tcl

# Paremeters:
set sys sin+dna
set conc 1.0
set psf ${sys}_hex.psf
set pdb ${sys}_hex.pdb
set outPrefix ${sys}_ions

# Load the system.
mol load psf $psf pdb $pdb

# Compute the number of ions to add.
set posSel [atomselect top "name POT"]
set posNum [$posSel num]
set negSel [atomselect top "name CLA"]
set negNum [$negSel num]
set sel [atomselect top "name OH2"]
set nw [$sel num]
set alpha 55.523
set ni [expr int(floor($conc*$nw/($alpha + 2.0*$conc)+0.5))]
$sel delete

# Get the charge.
set all [atomselect top all]
set other [atomselect top "not name POT CLA"]
set nq [expr int(floor([measure sumweights $other weight charge]+0.5))]
set nna [expr $ni - $posNum - $nq]
set ncl [expr $ni - $negNum]
puts "posNum0: $posNum"
puts "posNumQ: $nq"
puts "posNum1: $nna"
puts "negNum0: $negNum"
puts "negNum1: $ncl"
$all delete
$other delete
mol delete top

# Use autoionize.
package require autoionize
autoionize -psf $psf -pdb $pdb -nions [list "POT $nna" "CLA $ncl"] -o $outPrefix
