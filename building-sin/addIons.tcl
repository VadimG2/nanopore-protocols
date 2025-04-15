resetpsf
mol load psf sin_hex.psf pdb sin_hex.pdb
set conc 1.0
set water [atomselect top "name OH2"]
set num [expr {int(floor($conc*[$water num]/(55.523 + 2.0*$conc) + 0.5))}]
package require autoionize
autoionize -psf sin_hex.psf -pdb sin_hex.pdb -nions [list "POT $num" "CLA $num"] -o sin_ions
