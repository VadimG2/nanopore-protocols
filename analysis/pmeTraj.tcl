# Author: Jeff Comer <jcomer2@illinois.edu>
package require pmepot

proc compute {name structPrefix dcdList dcdFreq stride} {
    set displayPeriod 20
    set timestep 1.0
    set outDir .
    set outFile "$outDir/pot_${name}.dx"
    set pmeSize {96 96 192}

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc

    # Load the system.
    mol load psf $psf pdb $pdb
 
    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]
    }

    # Run pmePot.
    pmepot -mol top -xscfile $xsc -ewaldfactor 0.25 -grid $pmeSize -dxfile $outFile -frames all
    mol delete top
}
