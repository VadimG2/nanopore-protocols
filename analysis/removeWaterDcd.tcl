# Author: Jeff Comer <jcomer2@illinois.edu>

# Remove the directories from a path, leaving only the file name.
proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

# Remove water and silicon material from a dcd file.
proc removeWater {structPrefix dcd} {
    # Select the water or silicon material.
    set selText "(not water) and (not resname SIN SIO2)" 
    # Prefix to add to the output dcd files:
    set outPrefix "nw_"; 

    # Load the system.
    mol load psf $structPrefix.psf pdb $structPrefix.pdb
    set sel [atomselect top $selText]
    set dcdName [trimPath $dcd]
    set structName [trimPath $structPrefix]

    # Write the structure files for the resulting system.
    $sel writepsf nw_${structName}.psf
    $sel writepdb nw_${structName}.pdb

    # Load the trajectory.
    animate delete all
    mol addfile $dcd waitfor all
    set last [expr {[molinfo top get numframes]-1}]
    
    # Write the dcd files for the selection.
    animate write dcd "nw_$dcdName" beg 0 end $last waitfor all sel $sel top

    $sel delete
    mol delete top
}
