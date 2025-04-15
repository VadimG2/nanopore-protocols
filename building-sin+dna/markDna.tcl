set all [atomselect top all]
set sel [atomselect top "segname P1"]
$sel set beta 0.0
$sel set beta 1.0
$sel set occupancy 1.0
$all writepdb specific.pdb
