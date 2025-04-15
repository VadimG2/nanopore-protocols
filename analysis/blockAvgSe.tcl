#! /usr/bin/tclsh

# Compute block averages and standard deviations.
# Use with: tclsh blockAvgSe.tcl window_size input_size
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh blockAvgSe.tcl windowSize inputFile0 \[inputFile1\]..."
    exit
}

set window [lindex $argv 0]
# Input:
set inFileList [lrange $argv 1 end]
# Output:
set outSuffix $window.block

puts "Forming block averages with a window of size ${window}..."
if {![string is integer $window]} {
    puts "ERROR: The window size must be an integer."
    exit
}

foreach inFile $inFileList {
    # Open the files.
    set in [open $inFile r]
    set outFile $inFile.$outSuffix
    set out [open $outFile w]
    puts "Computing block averages for $inFile..."
    
    # Initialize the sums.
    set n 0
    set sumX 0.0
    set sumY 0.0
    set sumYY 0.0
    
    set nCol -1
    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} {continue}
	
	if {[string match "#*" $line]} {
	    puts $out $line
	} else {
	    set item [concat $line]
	    set nCol1 [llength $item]
	    if {$nCol < 0} {
		set nCol $nCol1
		puts "Reading $nCol columns."
		if {$nCol < 2} {
		    puts "ERROR: At least two columns are required."
		    exit
		}
		if {$nCol != 2} {puts "Warning: Ignoring [expr $nCol-2] columns."}
	    } elseif { $nCol1 != $nCol } {
		puts "Warning: Expecting $nCol data columns but found $nCol1!"
		continue
	    }
	    
	    # Add to the sums.
	    incr n
	    set x [lindex $item 0]
	    set sumX [expr {$sumX + $x}]
	    set y [lindex $item 1]
	    set sumY [expr {$sumY + $y}]
	    set sumYY [expr {$sumYY + $y*$y}]
	    
	    # Compute the averages for a completed window.
	    if {$n >= $window} {
		# Compute the statistics.
		set meanX [expr {$sumX/$n}]
		set meanY [expr {$sumY/$n}]
		set errY [expr {sqrt(($sumYY - $sumY*$sumY/$n)/($n-1)/$n)}] 

		# Write the results.
		puts $out "$meanX $meanY $errY"

		# Reinitialize the sums.
		set n 0
		set sumX 0.0
		set sumY 0.0
		set sumYY 0.0
	    }
	}
    }
    
    # Make a final data point if we have a lot of data left.
    if { $n > $window/2 + 1 } {
	puts "Making a final block average with the final $n data points..."
	set meanX [expr {$sumX/$n}]
	set meanY [expr {$sumY/$n}]
	set errY [expr {sqrt(($sumYY - $sumY*$sumY/$n)/($n-1)/$n)}] 

	# Write the results.
	puts $out "$meanX $meanY $errY"
    }

    close $out
    close $in
    puts "The block averages and standard errors were successfully written to `${outFile}'."
}



