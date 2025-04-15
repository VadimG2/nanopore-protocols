set sel [atomselect top "resname SIN"]
foreach quiet {0} { set pos [$sel get {x y z}] }
set out [open sin_positions.txt w]
foreach r $pos { puts $out $r }
close $out
