#!/usr/bin/tclsh
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh meanValueCut cutTime inputFile0 \[inputFile1\]..."
    exit
}

# Input:
set cutTime [lindex $argv 0]
set inFileList [lrange $argv 1 end]

puts "Analyzing data after t=$cutTime."
if {![string is double $cutTime]} {
    puts "ERROR: cutTime must be a real number."
    exit
}

set dimList {val}
foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    
    foreach dim $dimList {
	set sum($dim) 0
	set sumSq($dim) 0
    }
    set n 0

    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} {break}
	
	if {[string match "#*" $line]} {continue}
	set tok [concat $line]
	set x(tim) [expr double([lindex $tok 0])]
	set x(val) [expr double([lindex $tok 1])]

	if {$x(tim) > $cutTime} {
	    foreach dim $dimList {
		set sum($dim) [expr $sum($dim) + $x($dim)]
		set sumSq($dim) [expr $sumSq($dim) + $x($dim)*$x($dim)]
	    }
	    incr n
	}
    }
    
    puts "\nRead $inFile"
    foreach dim $dimList {
	set mean($dim) [expr $sum($dim)/$n]
	set err($dim) [expr sqrt(($sumSq($dim) - $sum($dim)*$sum($dim)/$n)/($n-1)/$n)]
	puts "mean $dim: $mean($dim) +/- $err($dim)"
    }
    
    close $in
}

