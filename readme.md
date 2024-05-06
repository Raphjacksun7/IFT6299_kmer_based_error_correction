
# K-mer Based Error Correction Tool

## Overview
This Python script corrects errors in DNA sequencing data by analyzing k-mer frequencies. Using Jellyfish for efficient k-mer counting, the tool dynamically adjusts corrections based on computed frequencies and quality scores. This approach optimizes nucleotide base selection throughout the correction process to enhance data quality for downstream analysis.


## Purpose
This tool aims to enhance the reliability of DNA sequencing data by reducing the frequency of erroneous k-mers. These corrections are crucial for improving the outcomes of downstream processes like genome assembly and SNP discovery.

## Algorithm
The script uses the "Weighted Adjusted Nucleotide Base Selection for" method:
1. **K-mer Counting**: Counts k-mers in the input FASTQ files using Jellyfish.
2. **Frequency Analysis**: Identifies k-mers with frequencies below a defined threshold (i0), suggesting potential errors.
3. **Error Correction**: For each potentially erroneous k-mer, the base with the lowest quality score is targeted for correction.
4. **Weighted Selection**: A dynamic weighting system is applied to select the most probable correct nucleotide based on existing data and past correction success.
5. **Adjustments**: Weights are adjusted based on the outcomes of corrections, employing an epsilon-greedy strategy to balance between exploring new corrections and exploiting known successful corrections.

## Requirements

- **Python 3.6 or higher**
- **Biopython**
- **NumPy**
- **Jellyfish** (must be accessible from the command line)

## Installation Guide

### Python Installation

Ensure Python 3.6 or later is installed on your system. Check your Python version with:

```bash
python3 --version
```

If Python is not installed, you can install it using your distribution's package manager. For Ubuntu, for example:

```bash
sudo apt update
sudo apt install python3 python3-pip
```

### Dependencies Installation

Install Biopython and NumPy using pip. These libraries are essential for handling biological data and performing numerical operations, respectively:

```bash
pip3 install biopython numpy
```

### Jellyfish Installation

Jellyfish is employed for rapid, memory-efficient k-mer counting in DNA sequences. Install it using your package manager:

```bash
sudo apt install jellyfish
```

## Usage

### Running the Script

Navigate to the directory containing the script and execute it with Python. Specify the necessary and optional arguments as detailed below:

```bash
python3 kmer_error_correction.py <input_fastq_file> -k <kmer_length> -d <output_directory> --evaluate --max-reads <max_number_of_reads>
```

### Command Line Arguments

- `<input_fastq_file>`: Path to the input FASTQ file.
- `-k <kmer_length>`: Specifies the length of k-mers (default is 17).
- `-d <output_directory>`: Designates the directory for storing output and results (default is `output`).
- `--evaluate`: Enables evaluation of correction efficiency (optional).
- `--max-reads <max_number_of_reads>`: Limits the number of reads processed, useful for quick testing or limited analysis (optional).

### Example Usage

```bash
python3 kmer_error_correction.py sample.fastq -k 21 -d results --evaluate --max-reads 5000
```

This command corrects errors in `sample.fastq` using a k-mer length of 21, saves results in the `results` directory, evaluates the correction efficiency, and processes only the first 5,000 reads.

## Output
- **Corrected FASTQ Files**: Saved in the specified output directory.
- **Evaluation Reports**: Generated if evaluation is enabled, comparing k-mer frequencies before and after corrections.

## License
This software is released under the [MIT License](LICENSE.md).

## Troubleshooting
Ensure Jellyfish is properly installed and accessible in your system's PATH. Refer to Python error logs for detailed troubleshooting guidance.

