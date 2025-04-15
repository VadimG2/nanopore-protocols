# Compute the mean value of the system length along z (c_z) for an NPT simulation.
# Scale the last frame of the equilibration to the mean value to produce a pdb.
# Produce an appropriate *.xsc file to go with it.
# Author: Jeff Comer <jcomer2@illinois.edu>

set sys sin+dna
set selText all
set minTime 5000; # ignore timesteps below this for average of c_z
# Input:
set structPrefix ../building-${sys}/${sys}_ions; # Your *.pdb, *.psf files
set inPrefix ../building-${sys}/${sys}_eq; # Your *.restart.coor and *.xst files
# Output:
set outPrefix scaled_${sys}_ions

mol delete all

proc writeXsc {xscList outFile} {
    set out [open $outFile w]
    puts $out "# NAMD extended system configuration restart file"
    puts $out "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
    puts $out $xscList
    close $out
    return
}

proc meanXst {inFile outPrefix minTime} {
    # Open the files.
    set in [open $inFile r]

    set valIndex 9
    set out [open $outPrefix.cz.dat w]
    set sum 0
    set sumSq 0
    set n 0

    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} {break}
	
	if {[string match "#*" $line]} {continue}
	set tok [concat $line]
	set tim [lindex $tok 0]
	set x [lindex $tok $valIndex]
	puts $out "$tim $x"

	if {$tim < $minTime} continue
	set sum [expr $sum + $x]
	set sumSq [expr $sumSq + $x*$x]
	
	incr n
    }

    puts "\nRead $inFile"
    set mean [expr $sum/$n]
    set err [expr sqrt(($sumSq - $sum*$sum/$n)/($n-1)/$n)]
    puts "mean size: $mean +/- $err"

    set finalX $x
    set xscList $tok
    lset xscList $valIndex $mean
    writeXsc $xscList $outPrefix.xsc

    close $out
    close $in

    return [list $mean $finalX]
}

# Calculate the mean system size.
set meanList [meanXst $inPrefix.xst $outPrefix $minTime]
set meanSize [lindex $meanList 0]
set lastSize [lindex $meanList 1]

mol load psf $structPrefix.psf
mol addfile $inPrefix.restart.coor
set scaleX 1.0
set scaleY 1.0
set scaleZ [expr {$meanSize/$lastSize}]

puts "\nNOTE:"
puts "original size: $lastSize"
puts "mean size: $meanSize"
puts "scale factor: $scaleZ"

# Scale the positions.
set sel [atomselect top $selText]
foreach zero {0} {set pos [$sel get {x y z}]}
set newPos {}
foreach r $pos {
    foreach {x y z} $r {break}
    set x [expr {$scaleX*$x}]
    set y [expr {$scaleY*$y}]
    set z [expr {$scaleZ*$z}]
    lappend newPos [list $x $y $z]
}
$sel set {x y z} $newPos
$sel delete

# Write the results.
set all [atomselect top all]
$all writepdb $outPrefix.pdb
$all delete
