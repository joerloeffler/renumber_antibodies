# Antibody PDB Renumbering Script

This script automates the process of renumbering antibody chains within a PDB file according to a standard numbering scheme (e.g., Kabat, Chothia). It uses ANARCI to identify and number the variable domains, renames the chains to 'H' (Heavy) and 'L' (Light), and merges the processed chains into a single, clean PDB file.

## Features

-   **Automatic Chain Identification**: Uses ANARCI to determine if a chain is an antibody and whether it is a heavy or light chain.
-   **Standardized Renumbering**: Renumbers residues in the variable domain according to a specified scheme (Kabat by default).
-   **Handles Non-Canonical Tails**: Correctly renumbers C-terminal residues that are not part of the canonical antibody domain by numbering them sequentially after the last ANARCI-numbered residue.
-   **Chain Renaming**: Renames the identified antibody chains to `H` and `L` for consistency.
-   **Clean Output**: Merges the renumbered chains into a single, new PDB file, leaving the original file untouched.

---

## Dependencies

To run this script, you will need the following software and Python libraries installed.

### Software
-   **Python** (Version 3.7 or higher recommended)
-   **ANARCI**: Must be installed and accessible from your system's PATH. You can find installation instructions [here](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/).

### Python Libraries
-   **BioPython**: A crucial library for handling PDB files and biological sequences.

---

## Installation

The recommended way to set up the environment is by using `conda` to ensure all dependencies are handled correctly.

1.  **Clone the repository:**
    ```bash
    git clone (https://github.com/joerloeffler/renumber_antibodies.git)
    cd renumber_antibodies
    ```

2.  **Create and activate a Conda environment:**
    This command creates a new environment named `renumAB` with Python 3.9.
    ```bash
    conda create -n renumAB python=3.9
    conda activate renumAB
    ```

3.  **Install BioPython:**
    ```bash
    pip install biopython
    ```

4.  **Verify ANARCI installation:**
    Make sure ANARCI is installed and working by running:
    ```bash
    ANARCI --help
    ```
    If this command fails, you need to install ANARCI and/or add it to your system's PATH.

---

## Usage

Run the script from your terminal, providing the path to your input PDB file.

#### Basic Usage
This will process the PDB file using the default **Kabat** numbering scheme.

```bash
python renumber.py /path/to/your/antibody.pdb
