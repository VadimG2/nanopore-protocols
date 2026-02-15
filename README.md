# Nanopore Simulation Automation Script (Single Peptide Chain)

This Bash script automates the setup and execution of molecular dynamics simulations for a **single peptide chain** translocating through a **Si₃N₄ nanopore** under an electric field. It is based on the tutorial *“Modeling Nanopores for Sequencing DNA”* but adapted for a peptide. The script handles:

- Building the peptide structure (Section 3.1)
- Building the pure Si₃N₄ nanopore system (Section 3.5)
- Simulating both the empty pore and the pore+peptide system under bias (Section 3.9)
- (Optionally) building the combined pore+peptide system (Section 3.6) – currently commented out but can be enabled.

All steps are executed in the appropriate directories, and VMD Tcl scripts are run via stdin piping to ensure clean exits.

---

## Directory Structure

The script expects the following layout inside your working directory:

```
.
├── building-dna/          # Peptide building (Section 3.1)
│   ├── separate.tcl
│   ├── dsdna.pdb          (initial peptide structure, expected by separate.tcl)
│   └── ...
├── building-sin/          # Pure Si₃N4 pore building (Section 3.5)
│   ├── sin.psf / sin.pdb  (initial pore files)
│   ├── solvate.tcl
│   ├── cutWaterHex.tcl
│   ├── addIons.tcl
│   ├── markRestraints.tcl
│   ├── sin_min.namd
│   ├── sin_eq.namd
│   └── ...
├── building-sin+dna/      # Combined pore+peptide building (Section 3.6) – commented out by default
│   ├── combine.tcl
│   ├── adjustPos.tcl
│   ├── defineRestraints.tcl
│   ├── writePos.tcl
│   ├── markDna.tcl
│   ├── sin+dna_min.namd
│   ├── sin+dna_eq.namd
│   └── ...
├── running-sin/           # Simulation runs (Section 3.9)
│   ├── scaleToMeanNptSize.tcl
│   ├── sin_20V.namd
│   ├── sin+dna_20V.namd   (commented out)
│   └── ...
├── c32b1/toppar/          # Topology files (path expected by autopsf step)
│   └── top_all27_prot_na.rtf
└── grid/                  # External force grid utility (optional, for peptide-specific forces)
    └── thirdForce
```

You must create these directories and populate them with the required input files (PSF/PDB for the pore, Tcl scripts, NAMD configuration files) before running the script.

---

## Prerequisites

- **VMD** – for structure manipulation and system building.
- **NAMD3** – for MD simulations (tested with `namd3`).
- **Topology files** – e.g., `top_all27_prot_na.rtf` (CHARMM27) placed in `c32b1/toppar/`.
- **Bash** – the script uses `set -e` and standard Unix tools.

---

## Configuration

Before running, you **must** customise a few settings:

### 1. Ion Concentration
Edit `addIons.tcl` in the respective building directory (`building-sin` or `building-sin+dna`). Look for a line like:
```tcl
set conc 1.0
```
Change the value to your desired salt concentration (e.g., 0.15 M).

### 2. Peptide Charge
The script assumes a **negatively charged peptide**. If your peptide is positively charged, you **must** modify the solvation command in **Section 3.6** (if you enable it) and possibly in **Section 3.5**.  
Look for the comment block at the top of the script:
```bash
###########>>> WARNING! THIS SCRIPT FOR NEGATIVELY CHARGED PEPTIDES! IF YOU WANT POSITIVE CHARGED TO BE SIMULATED, CHANGE: <<<##############################
```
Follow the instructions there to swap the `-minmax` values in the `solvate` command.

### 3. Simulation Box Dimensions
The `-minmax` arguments in solvate commands (e.g., `{{-55 -55 -97} {55 55 167}}`) must match your system’s dimensions. Adjust them according to your pore and peptide size.

### 4. Scale Script Mode
The script `scaleToMeanNptSize.tcl` in `running-sin` is edited **in-place** to switch between `sin` (empty pore) and `sin+dna` (pore+peptide) modes. The script currently uses a pattern substitution; ensure the patterns match those in your Tcl file.

---

## Usage

1. Place all required input files in the directories as described.
2. Make the script executable:
   ```bash
   chmod +x automate.sh
   ```
3. Run it from the **base directory** (the one containing `building-dna`, `building-sin`, etc.):
   ```bash
   ./automate.sh
   ```

The script will:
- Build the peptide structure (`building-dna`)
- Build the pure Si₃N₄ pore system (`building-sin`)
- Run NAMD simulations for the empty pore and (if uncommented) the pore+peptide system under a 20 V bias (`running-sin`)

---

## What the Script Does (Step by Step)

### Section 3.1 – Peptide Building
- Runs `separate.tcl` (expected to produce `dsdna_A.pdb` and `dsdna_B.pdb`).
- Renames `dsdna_A.pdb` → `dsdna.pdb` and deletes the empty `dsdna_B.pdb`.
- Uses `autopsf` to generate a PSF file for the peptide (`dsdna.psf`) with the CHARMM27 topology.

### Section 3.5 – Pure Si₃N₄ Pore Building
- Solvates the pore (`sin.psf`/`sin.pdb`) using `solvate`.
- Trims water to periodic boundaries (`cutWaterHex.tcl`).
- Adds ions (`addIons.tcl`).
- Defines harmonic restraints for Si₃N₄ atoms (`markRestraints.tcl`).
- Sets Langevin thermostat atoms (`beta=1` for Si₃N₄).
- Runs a short minimization (`sin_min.namd`) and NPT equilibration (`sin_eq.namd`).

### Section 3.6 – Pore+Peptide Building (Commented Out by Default)
If you uncomment this section, it will:
- Combine the pore and peptide (`combine.tcl`).
- Adjust the peptide position (`adjustPos.tcl`).
- Solvate, cut water, add ions, define restraints, and set up thermostat atoms.
- Generate a force grid for the peptide using an external utility (`../grid/thirdForce`).
- Run minimization and NPT equilibration.

### Section 3.9 – Simulations under Electric Field
- Switches to the `running-sin` directory.
- Runs `scaleToMeanNptSize.tcl` to scale the equilibrated system to the target size (this script is expected to read restart files from the building step and produce scaled coordinates).
- Launches NAMD production runs:
  - `sin_20V.namd` for the empty pore.
  - (Commented) `sin+dna_20V.namd` for the pore+peptide system.

---

## Important Notes

- The script uses `set -e` – it will **stop immediately** if any command fails. This helps catch errors early.
- VMD scripts are run by piping their content into `vmd -dispdev text`; ensure that each Tcl script ends gracefully (no interactive prompts).
- Some internal errors in `scaleToMeanNptSize.tcl` (e.g., division by zero, missing XSC files) are **not** handled by the script – you must check the logs.
- If you enable Section 3.6, verify that the `thirdForce` executable exists and that `sin_positions.txt` is generated correctly.
- The NAMD commands use `+p16 +devices 0` (adjust to your system). Modify `NAMD_OPTS` in the script to match your available CPUs/GPUs.

---

## Troubleshooting

- **Pattern not found in scale script**: If editing `scaleToMeanNptSize.tcl` fails, check that the patterns (`sin` / `sin+dna`) exactly match those in the Tcl file.
- **Missing restart files after equilibration**: Ensure the NAMD config files write `.restart.coor`, `.restart.vel`, and `.restart.xsc`.
- **autopsf cannot find topology**: Verify the path `../c32b1/toppar/top_all27_prot_na.rtf` exists relative to `building-dna`.

---

## License

This script is provided as-is for academic use. Modify and distribute freely, but please cite the original tutorial and the authors of VMD/NAMD.
