#!/bin/bash

# === Automation Script for Nanopore Simulation Protocol ===
# Based on "Modeling Nanopores for Sequencing DNA" Tutorial PDF
# Modified for SINGLE PEPTIDE CHAIN simulation.
# Automates Sections: 3.1 (Peptide build), 3.6 (Combine Si3N4+Peptide), 3.9 (Run Si3N4 Systems)
# Uses autopsf for Peptide PSF generation.
# Runs VMD Tcl scripts via stdin piping to ensure exit.
# Edits scaleToMeanNptSize.tcl in place to switch between 'sin' and 'sin+dna' modes.
# Executes scaling and NAMD simulation steps strictly from the 'running-sin' directory.
# Removes argument passing for paths to scaleToMeanNptSize.tcl, assuming it handles paths internally.
# NOTE: Internal errors in scaleToMeanNptSize.tcl (divide by zero, missing xsc write) likely remain.

############################################################################################################################################################
###########>>> WARNING! THIS SCRIPT FOR NEGATIVELY CHARGED PEPTIDES! IF YOU WANT POSITIVE CHARGED TO BE SIMULATED, CHANGE: <<<##############################

# package require solvate; solvate sin+dna.psf sin+dna_placed.pdb -minmax {{-55 -55 -97} {55 55 167}} -o sin+dna_sol; exit <<<##############################

# to THIS:

# package require solvate; solvate sin+dna.psf sin+dna_placed.pdb -minmax {{-55 -55 -197} {55 55 97}} -o sin+dna_sol; exit <<<#############################

############################################################################################################################################################

###########>>> ATTENTION! HOW TO CHANGE ION CONCENTRATIONS <<<##############################################################################################

# You must proceed to building-sin+dna directory and change addIons.tcl script. You must change "set conc 1.0" to you preferable concentration!

#############################################################################################################################################################

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---
VMD_EXEC="vmd"
NAMD_EXEC="namd3"
NAMD_OPTS="+p16 +devices 0"
BASE_DIR=$(pwd)
SCALE_SCRIPT="scaleToMeanNptSize.tcl" # Script name relative to running-sin

# Function to run VMD script via stdin piping (No Tcl args expected)
run_vmd_script() {
  local script_file="$1"
  if [ ! -f "$script_file" ]; then echo "ERROR: VMD script not found: $script_file in $(pwd)"; exit 1; fi
  echo "Running VMD script: $script_file (via stdin)"
  # Pipe script content and exit command to VMD
  (cat "$script_file"; echo "exit") | ${VMD_EXEC} -dispdev text
  echo "Finished VMD script: $script_file"
}

# Function to safely edit the scale script
edit_scale_script() {
    local from_pattern="$1"; local to_pattern="$2"; local script_path="$3"
    echo "Editing $script_path: Changing '$from_pattern' to '$to_pattern'"
    if ! grep -q "$from_pattern" "$script_path"; then echo "WARNING: Pattern '$from_pattern' not found in $script_path. Edit may have already occurred or pattern differs."; return 1; fi
    # Use ::: as delimiter in sed to avoid conflicts with paths containing /
    sed "s:::$from_pattern:::$to_pattern::" "$script_path" > "${script_path}.tmp" && mv "${script_path}.tmp" "$script_path"
    if [ $? -ne 0 ]; then echo "ERROR: Failed to edit $script_path using sed."; exit 1; fi
    if ! grep -q "$to_pattern" "$script_path"; then echo "ERROR: Verification failed. Pattern '$to_pattern' not found after editing $script_path."; exit 1; fi
    echo "Successfully edited $script_path."; return 0
}

echo "=== Starting Nanopore Simulation Automation Script (PEPTIDE MOD + Script Edit + Dir Fix V2) ==="
# ... (Остальная информация) ...

# === Section 3.1: Building Peptide Structure ===
echo "--- Section 3.1: Building Peptide Structure ---"
cd "${BASE_DIR}/building-dna"; echo "[3.1] Current directory: $(pwd)"
echo "[3.1 Step 3 - Modified] Running separate.tcl ..."; run_vmd_script separate.tcl
echo "[3.1 Post-Step 3] Renaming peptide PDB and cleaning up..."
if [ -f dsdna_A.pdb ]; then mv dsdna_A.pdb dsdna.pdb; echo "Renamed dsdna_A.pdb to dsdna.pdb"; else echo "ERROR: dsdna_A.pdb not found"; exit 1; fi
if [ -f dsdna_B.pdb ]; then rm dsdna_B.pdb; echo "Removed empty dsdna_B.pdb"; fi
echo "[3.1 Step 4 - Modified] Generating Peptide PSF file using autopsf..."
cat << EOF > run_autopsf_peptide.tcl
package require autopsf; package require topotools
mol load pdb dsdna.pdb
set topology_dir "../c32b1/toppar"; set topology_file "top_all27_prot_na.rtf"; set topology_path "\$topology_dir/\$topology_file"
if {![file exists \$topology_path]} { puts stderr "ERROR: Topology file not found at \$topology_path"; exit 1 }; set topolfiles [list \$topology_path]
autopsf -top \$topolfiles
set sel [atomselect top all]; \$sel writepsf dsdna.psf; \$sel writepdb dsdna.pdb
if {![file exists dsdna.psf]} { puts stderr "ERROR: autopsf failed to create dsdna.psf"; exit 1 }
puts "autopsf completed successfully for peptide using \$topology_file."
exit
EOF
${VMD_EXEC} -dispdev text -e run_autopsf_peptide.tcl; rm run_autopsf_peptide.tcl
echo "[3.1 Step 4] dsdna.psf and updated dsdna.pdb generated."
cd "${BASE_DIR}"; echo "--- Section 3.1 Finished ---"; echo ""

# === Section 3.6: Building the synthetic pore-PEPTIDE system ===
echo "--- Section 3.6: Building Synthetic Pore-Peptide System ---"
cd "${BASE_DIR}/building-sin+dna"; echo "[3.6] Current directory: $(pwd)"
echo "[3.6 Step 1] Combining Si3N4 nanopore and PEPTIDE..."; run_vmd_script combine.tcl
echo "[3.6 Step 2] Adjusting PEPTIDE position..."; run_vmd_script adjustPos.tcl
echo "[3.6 Step 3] Solvating the combined system..."; cat << EOF > run_vmd_solvate_sin.tcl
package require solvate; solvate sin+dna.psf sin+dna_placed.pdb -minmax {{-55 -55 -97} {55 55 167}} -o sin+dna_sol; exit
EOF
${VMD_EXEC} -dispdev text -e run_vmd_solvate_sin.tcl; rm run_vmd_solvate_sin.tcl
echo "[3.6 Step 4] Cutting water to periodic boundaries..."; run_vmd_script cutWaterHex.tcl
echo "[3.6 Step 5] Adding ions..."; run_vmd_script addIons.tcl
echo "[3.6 Step 6] Defining harmonic restraints for Si3N4..."; run_vmd_script defineRestraints.tcl
mv sin+dna_restrain.pdbexit sin+dna_restrain.pdb
echo "[3.6 Step 7] Setting up membrane thermostat atoms..."; cat << EOF > run_thermostat_setup.tcl
mol load psf sin+dna_ions.psf pdb sin+dna_ions.pdb; set all [atomselect top all]; set sel [atomselect top "resname SIN"]; \$all set beta 0.0; \$sel set beta 1.0; \$all writepdb sin+dna_langevin.pdb; exit
EOF
${VMD_EXEC} -dispdev text -e run_thermostat_setup.tcl; rm run_thermostat_setup.tcl
echo "[3.6 Step 8] Writing membrane atom positions..."; cat << EOF > run_writePos_wrapper.tcl
mol load psf sin+dna_ions.psf pdb sin+dna_ions.pdb; source writePos.tcl; exit
EOF
echo "Running writePos.tcl via wrapper..."; ${VMD_EXEC} -dispdev text -e run_writePos_wrapper.tcl; rm run_writePos_wrapper.tcl; echo "Finished writePos.tcl (via wrapper)"
echo "[3.6 Step 9] Marking PEPTIDE atoms..."; cat << EOF > run_markDna_wrapper.tcl
mol load psf sin+dna_ions.psf pdb sin+dna_ions.pdb; source markDna.tcl; exit
EOF
echo "Running markDna.tcl via wrapper..."; ${VMD_EXEC} -dispdev text -e run_markDna_wrapper.tcl; rm run_markDna_wrapper.tcl; echo "Finished markDna.tcl (via wrapper)"
echo "[3.6 Step 10] Generating PEPTIDE-specific force grid..."; if [ -f ../grid/thirdForce ]; then if [ -s sin_positions.txt ]; then ../grid/thirdForce sin_positions.txt grid_basis.txt 1 2 2 specific2-2.dx; else echo "ERROR: sin_positions.txt empty"; exit 1; fi else echo "WARNING: ../grid/thirdForce not found"; fi
echo "[3.6 Step 11] Minimizing Si3N4+Peptide system..."; ${NAMD_EXEC} ${NAMD_OPTS} sin+dna_min.namd > sin+dna_min.log
echo "[3.6 Step 12] Equilibrating Si3N4+Peptide system (NPT)..."; ${NAMD_EXEC} ${NAMD_OPTS} sin+dna_eq.namd > sin+dna_eq.log
cd "${BASE_DIR}"; echo "--- Section 3.6 Finished ---"; echo ""


# === Section 3.9: Simulating synthetic nanopore with PEPTIDE under E-field ===
echo "--- Section 3.9: Simulating Si3N4 Systems with E-Field ---"
# !!! ПЕРЕХОДИМ В РАБОЧУЮ ДИРЕКТОРИЮ ДЛЯ ЭТОЙ СЕКЦИИ !!!
cd "${BASE_DIR}/running-sin"
echo "[3.9] Current directory: $(pwd)"

# Запускаем VMD скрипт (он должен сам найти пути типа ../building-sin/...)
#run_vmd_script scaleToMeanNptSizeSIN.tcl # Аргументы не передаем

# --- Step 4a: Run NAMD Simulation for Empty Pore ---
#echo "[3.9 Step 4a] Simulating EMPTY PORE under 20V bias (NVT)..."
# Убедитесь, что sin_20V.namd читает scaled_sin_ions.pdb/xsc/etc из ТЕКУЩЕЙ директории
#namd3 +devices 0 +p16 sin_20V.namd > sin_20V.log
#echo "Running: $NAMD_CMD_EMPTY"

#Step 3 (Part 2): Scale Pore+Peptide System ---
echo "[3.9 Step 3b] Scaling system size for PORE+PEPTIDE (using 'set sys sin+dna')..."
# Запускаем VMD скрипт (он должен сам найти пути типа ../building-sin+dna/...)
run_vmd_script scaleToMeanNptSizeSINDNA.tcl # Аргументы не передаем
# Tcl скрипт должен создать scaled_sin+dna_ions.pdb/xsc в текущей директории (running-sin)

# --- Step 4b: Run NAMD Simulation for Pore+Peptide System ---
echo "[3.9 Step 4b] Simulating PORE+PEPTIDE system under 20V bias (NVT)..."
# Убедитесь, что sin+dna_20V.namd читает scaled_sin+dna_ions.pdb/xsc/etc из ТЕКУЩЕЙ директории
namd3 +devices 0 +p2 sin+dna_20V.namd > sin+dna_20V.log
echo "Running: $NAMD_CMD_PEPTIDE"

cd "${BASE_DIR}"
echo "--- Section 3.9 Finished ---"
echo ""

echo "=== Automation Script Finished ===="
echo "NOTE: If NAMD steps were skipped, please check the ERROR messages above regarding ${SCALE_SCRIPT}."

exit 0
