# seqFISH - Sequential RNA FISH Analysis Pipeline

This repository contains tools for designing and analyzing sequential RNA FISH experiments. It includes codebook generation, spot detection, and gene identification pipelines.


## Prerequisites

Before you start, make sure you have the following installed:

### 1. Git
Download and install Git from: https://git-scm.com/download

**Why Git?** Git is a version control system that allows you to download this entire repository with a single command.

### 2. Miniconda (Python Package Manager)
Download and install Miniconda from: https://www.anaconda.com/docs/getting-started/miniconda/install

**Why Miniconda?** It makes it easy to manage Python packages and create isolated environments for this project.

### 3. Visual Studio Code (Recommended)
Download and install VSCode from: https://code.visualstudio.com/

**Why VSCode?** It's a lightweight, beginner-friendly code editor with excellent Python support.

Once installed, install the **Python Extension**:
- Open VSCode
- Go to Extensions (Ctrl+Shift+X on Windows, Cmd+Shift+X on Mac)
- Search for "Python" and click "Install" on the Microsoft Python extension


## Quick Start

### Step 1: Open the Miniconda Terminal

**IMPORTANT:** The following commands ONLY work in the Miniconda/Anaconda terminal, NOT in regular Command Prompt or PowerShell!

**Windows:**
- Open the Start Menu
- Search for "Anaconda Prompt" or "Miniconda Prompt"
- Click to open it

**Or in VSCode:**
- Press Ctrl+` (backtick) to open the terminal
- At the bottom, click the dropdown menu and select "Command Prompt" or "Anaconda Prompt"
- Or manually select it: Terminal → New Terminal, then change the shell

**macOS/Linux:**
- Open the Terminal application

**In VSCode:**
- Press Ctrl+` (backtick) or go to Terminal → New Terminal

### Step 2: Clone the Repository

Copy and paste this command into your **Miniconda terminal** and press Enter:

```bash
git clone https://github.com/HafnerLab/seqFISH.git
cd seqFISH
```

This downloads the entire seqFISH repository to your computer.

### Step 3: Create a Conda Environment

Conda environments are like isolated "containers" for Python projects. This keeps dependencies organized.

Copy and paste this command into your **Miniconda terminal** and press Enter:

```bash
conda create -n seqfish_hafner -c conda-forge python=3.12
```

**This command MUST be run in the Miniconda/Anaconda terminal!**

Press `y` and Enter when asked to confirm.

This may take 1-2 minutes.

### Step 4: Activate the Environment

Every time you want to use seqFISH, activate this environment first by copying and pasting this into your **Miniconda terminal**:

```bash
conda activate seqfish_hafner
```

You should see `(seqfish_hafner)` at the start of your terminal prompt.

### Step 5: Install Required Packages

Copy and paste this command into your **Miniconda terminal** and press Enter:

```bash
pip install -r seqFISH/requirements.txt
```

This installs all necessary Python packages. It may take 2-5 minutes.

### You're Ready!

You have successfully set up seqFISH!

**Note:** Next time you want to use seqFISH, just:
1. Open the Miniconda terminal
2. Navigate to your seqFISH folder: `cd path/to/seqFISH`
3. Activate the environment: `conda activate seqfish_hafner`
4. Run your script


## Usage Guide

Each Python script performs a specific step in the analysis pipeline. Run them in order:

**First, make sure your Miniconda terminal is open and the environment is activated:**
```bash
conda activate seqfish_hafner
```

### A. Generate Readout Probe Sequences (Optional)

Generate candidate readout sequences for BLAST screening:

```bash
python A_Generate_readout_sequences_forBLAST.py
```

**When to use:** When designing new probe sets for your genes.

**Configuration in script:**
- `SEQ_LENGTH`: Length of readout probes (default: 20 nt)
- `NUM_SEQ`: Number of sequences to generate (default: 100)
- `MAX_OVERLAP`: Maximum overlap allowed (default: 5 nt)

---

### B. Filter BLAST Results (Optional)

Filter sequences based on specificity and overlap:

```bash
python B_Filter_BLAST_results.py
```

**When to use:** After running BLAST on your candidate sequences.

**Input required:** BLAST XML output and FASTA file

---

### C. Create Codebook

Generate hybridization codebook for your genes:

**Option 1: Command Line (Advanced)**
```bash
python helper_functions.py
```

**Option 2: Web Interface (Recommended for Beginners)**
```bash
python Codebook_Generator.py
```
Then open your web browser and go to: http://localhost:5000

**Supported modes:**
- **MERFISH**: Sequential FISH with readout probes
- **RAEFISH**: Rolling circle amplification FISH

---

### D. Channel Registration (Optional)

Perform channel-to-channel registration (aligns different color channels):

```bash
jupyter notebook D_NanoPyx_ChannelRegistration.ipynb
```

A Jupyter notebook window will open in your browser. Just follow the on-screen instructions.

**When to use:** If your microscope images have misaligned color channels.

---

### E. Drift Correction (Optional)

Correct for sample movement during imaging:

```bash
jupyter notebook E_NanoPyx_DriftCorrection.ipynb
```

**When to use:** If your sample moved during time-lapse imaging.

---

### F. Spot Detection and Gene Identification (Main Analysis)

The main analysis script. Detects RNA spots and identifies genes:

```bash
python F_FISH_Quant_detection.py
```

**Features:**
- Multi-format image support (TIFF, ND2, CZI, LIF)
- Automatic image dimension detection
- Edge exclusion and stitching artifact removal
- Interactive threshold adjustment
- Gene identification from codebook
- Quality control repeat validation
- High-resolution visualization

**What happens:**
1. A dialog appears asking you to select your image file(s)
2. You select your codebook Excel file
3. The script detects spots and identifies genes
4. Results are saved to an `output_[date]_[time]` folder

---

## Project Structure

```
seqFISH/
├── doc/
│   └── MOUSE_sequences.xlsx           # Mouse gene sequences and primers
├── A_Generate_readout_sequences_forBLAST.py
├── B_Filter_BLAST_results.py
├── helper_functions.py
├── D_NanoPyx_ChannelRegistration.ipynb
├── E_NanoPyx_DriftCorrection.ipynb
├── F_FISH_Quant_detection.py          # Main analysis script
├── Codebook_Generator.py              # Web interface
├── requirements.txt
├── README.md
├── templates/                         # Web interface templates
│   └── index.html
└── .gitignore
```

---

## System Requirements

- **Python**: 3.12 or newer (installed via Miniconda)
- **Memory**: 8 GB RAM minimum (16 GB recommended for large images)
- **Disk Space**: 10 GB minimum
- **Operating System**: Windows, macOS, or Linux

---

## Output Files

After running `F_FISH_Quant_detection.py`, these files are created in the `output_` folder:

| File | Description |
|------|-------------|
| `{prefix}_colocalized.csv` | Raw spot groups detected |
| `{prefix}_colocalized_identified.csv` | Spots with identified genes (main results) |
| `{prefix}_spots_indexing.csv` | Spot location reference |
| `{prefix}_seqFISH_summary.csv` | Summary statistics |
| `{prefix}_gene_colored_scatter.png` | Visualization (open in any image viewer) |
| `{prefix}_thresholds.csv` | Detection thresholds (can be reused for similar images) |

---

## Troubleshooting

### "Command not found: conda"
**Problem:** You're not using the Miniconda/Anaconda terminal.

**Solution:**
1. Close your current terminal
2. Open the "Anaconda Prompt" or "Miniconda Prompt" application (search in Start Menu)
3. Try your command again

---

### "Command not found: git"
**Problem:** Git is not installed or not in your system PATH.

**Solution:**
1. Download Git from https://git-scm.com/download
2. Install it (use default options)
3. Restart your Miniconda terminal

---

### "aicsimageio not installed"
**Problem:** The image reader package failed to install.

**Solution:**
```bash
conda activate seqfish_hafner
pip install aicsimageio[nd2,czi]
```

---

### "Cannot read ND2/CZI/LIF files"
**Problem:** Your specific image format needs additional support.

**Solution:**
```bash
conda activate seqfish_hafner
pip install nd2
pip install python-bioformats
```

---

### Tkinter dialogs not appearing (macOS)
**Problem:** Tkinter (for dialog windows) is missing.

**Solution:**
```bash
conda activate seqfish_hafner
conda install tk
```

---

### "MOUSE_sequences.xlsx not found"
**Problem:** The mouse sequence file is missing.

**Solution:**
Make sure `MOUSE_sequences.xlsx` is in the `doc/` folder.

---

### Environment activation not working
**Problem:** `conda activate` doesn't work in your terminal.

**Solution:**
1. Make sure you're using "Anaconda Prompt" or "Miniconda Prompt" (NOT regular Command Prompt)
2. On Windows, you may need to use: `activate seqfish_hafner` (without the `conda` prefix)
3. If still not working, reinstall Miniconda and check "Add to PATH" during installation

---

## References

- **Big-FISH**: Imbert et al. (2022) - https://github.com/fish-quant/big-fish
- **NanoPyx**: Laine et al. (2019) - https://github.com/HenriquesLab/NanoPyx


---

## Contact & Support

- **Issues/Questions:** https://github.com/HafnerLab/seqFISH/issues
- **Email:** wouter.droogers@donders.ru.nl

---

## Acknowledgments

This project builds on excellent work from:
- **Mueller Lab** - Big-FISH spot detection
- **Henriques Lab** - NanoPyx image processing
- **Hafner Lab** - Integration and workflow
