# seqFISH Codebook Generator

A web-based tool for generating hybridization codebooks for sequential RNA FISH experiments (MERFISH and RAEFISH).

## Prerequisites

- [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install)

## Installation

**1. Download this repository**

Click the green **Code** button on GitHub and select **Download ZIP**. Extract the folder to a location of your choice.

Or clone with Git (requires [Git](https://git-scm.com/download)):
```bash
git clone https://github.com/HafnerLab/seqFISH.git
```

Or download and extract via the command line without Git:

*Windows (PowerShell):*
```powershell
Invoke-WebRequest -Uri https://github.com/HafnerLab/seqFISH/archive/refs/heads/main.zip -OutFile seqFISH.zip
Expand-Archive seqFISH.zip -DestinationPath .
```

*macOS / Linux:*
```bash
curl -L https://github.com/HafnerLab/seqFISH/archive/refs/heads/main.zip -o seqFISH.zip
unzip seqFISH.zip
```

**2. Open the Miniconda terminal**

On Windows, open the **Anaconda Prompt** from the Start Menu. Navigate to the extracted folder:
```bash
cd path\to\iST-Codebook-Generator
```

**3. Create and activate a Python environment**

```bash
conda create -n codebook_generator -c conda-forge python=3.12
conda activate codebook_generator
```

**4. Install required packages**

```bash
pip install -r requirements.txt
```

## Running the application

With the environment activated, start the web interface:

```bash
conda activate codebook_generator
python Codebook_Generator.py
```

Then open your browser and go to: http://localhost:5000

## Using the application

### Mode

Select the experiment type:

- **MERFISH** — Sequential FISH using readout probes. Primer and readout sequences are loaded automatically from `doc/MOUSE_sequences.xlsx`.
- **RAEFISH** — Rolling circle amplification FISH. Adapter, readout probe, and gene sequences are loaded from `doc/MOUSE_sequences.xlsx`.

### Barcode type

- **Combinatorial** — Each gene is assigned a unique combination of readout numbers across hybridization rounds. This is the standard approach.
- **Binary** — Each gene is assigned a binary bit code (0/1) per imaging round, distributed across channels. Useful for designs where each channel is imaged independently per round.

### Gene input

- **Gene list file (Excel)** — Upload an `.xlsx` file with gene names in the first column. A second column with expression values is optional; if provided, the codebook will be balanced by expression level.
- **Number of genes (MERFISH only)** — Enter a count to generate a codebook with synthetic gene names, without requiring a gene list file.

RAEFISH always requires a gene list file, since sequences are looked up per gene from `MOUSE_sequences.xlsx`.

### Encoding probe processing (MERFISH only)

- **Reverse complement** — Applies reverse complement to the readout sequences in the encoding probes. Use this for **ProbeDealer**.
- **Leave intact** — Keeps the readout sequences as-is. Use this for **PaintSHOP**.

### Imaging rounds and channels

Set the number of hybridization/imaging rounds and fluorescence channels. The total number of addressable readouts equals rounds × channels in case of a combinatorial codebook.

### Balanced codebook

- **Enabled (default)** — Optimizes the codebook so each readout appears approximately equally often across genes, minimizing optical crowding. If expression values are provided, balancing is weighted by expression level.
- **Disabled** — Maximizes the number of valid barcodes without frequency optimization.

### Enforce pairwise uniqueness

- **Enabled (default)** — Strict mode: for any two barcodes, the union of their active readout positions must be unique. This reduces the risk of barcode confusion during decoding.
- **Disabled** — Relaxed mode: only requires a Hamming distance of at least 2 between any two barcodes.

### Output

Enter a file name for the output (e.g., `Codebook_output.xlsx`) and select the folder where it should be saved. Click **Generate Codebook** to start. Generation runs in the background; a status message will appear when complete.

## Output file

The generated Excel file contains multiple sheets depending on the mode:

| Sheet | Description |
|-------|-------------|
| `Codebook` | Gene-to-readout assignments |
| `Pattern` | Readout/adapter sequences per hybridization round and channel |
| `Secondary_sequences` | Encoding probe sequences (5′ and 3′) |
| `Ordering` | Padlock and priming sequences (RAEFISH only) |

## Project structure

```
├── Codebook_Generator.py      # Web interface (Flask)
├── helper_functions.py        # Codebook generation logic
├── templates/index.html       # Web interface template
├── doc/MOUSE_sequences.xlsx   # Primer, readout, and gene sequences
└── requirements.txt
```

## Troubleshooting

### "MOUSE_sequences.xlsx not found"
Make sure `MOUSE_sequences.xlsx` is present in the `doc/` folder.

### Tkinter folder dialog not appearing (macOS)
```bash
conda activate codebook_generator
conda install tk
```

### "conda activate" not working
Make sure you are using the **Anaconda Prompt**, not a regular terminal. On Windows, search for "Anaconda Prompt" in the Start Menu.

---

## Contact

- **Issues:** https://github.com/HafnerLab/iST-Codebook-Generator/issues
- **Email:** wouter.droogers@donders.ru.nl
