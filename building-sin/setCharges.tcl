# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
# Set the charges of silicon and nitrogen based
# to obtain a neutral system.

set silText "type SI"
set nitText "type NSI"
set silCharge 0.7679
set nitCharge -0.575925
set selTextLeaving "type SI and numbonds == 4"
#Input:
set psf sin_pore_types.psf
set pdb sin_pore_types.pdb
#Output:
set finalPsf sin_pore_charges.psf
set finalPdb sin_pore_charges.pdb

# Load the molecule.
mol load psf $psf pdb $pdb
set all [atomselect top all]

# Set the charges.
set sil [atomselect top $silText]
$sil set charge $silCharge
$sil delete
set sel [atomselect top $nitText]
$sel set charge $nitCharge

# Get the initial charge.
set charge [measure sumweights $all weight charge]
puts "Initial charge: $charge"

# Distribute it among selText.
set qi [lindex [$sel get charge] 0]
set num [$sel num]
set q [expr -$charge/$num + $qi]
set err [expr ($q-$qi)/$qi]
puts "Shifting the charge of `$nitText' from $qi to $q"
puts "Relative error: $err"
$sel set charge $q

# Write the intermediate result.
$all writepsf $finalPsf

# Clean up.
$sel delete
$all delete
mol delete top

# Reload to fix round off error.
mol load psf $finalPsf pdb $pdb
set all [atomselect top all]
set newCharge [measure sumweights $all weight charge]
puts "New charge: $newCharge"
set sel [atomselect top $selTextLeaving]
set q [expr [lindex [$sel get charge] 0] - $newCharge]
set index [lindex [$sel get index] 0]
$sel delete

# Shift.
puts "Shifting charge of atom $index by [expr -$newCharge]."
set sel [atomselect top "index $index"]
$sel set charge $q
set newCharge [measure sumweights $all weight charge]
puts "Final charge: $newCharge"

# Write the results.
$all writepdb $finalPdb
$all writepsf $finalPsf

$sel delete
$all delete
mol delete top
