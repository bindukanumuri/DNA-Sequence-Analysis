# DNA Sequence Analysis using BioPython

A Python-based bioinformatics project that performs DNA sequence analysis using the BioPython library. The project reads DNA sequences from FASTA files, analyzes nucleotide composition, calculates GC and AT content, performs transcription and translation, detects Open Reading Frames (ORFs), identifies restriction enzyme sites, and compares multiple DNA sequences using graphical visualization.

---

## Features

- Read DNA sequences from FASTA files
- Analyze single and multiple DNA sequences
- Calculate sequence length
- Count nucleotide frequencies (A, T, G, C)
- Calculate GC Content
- Calculate AT Content
- DNA to RNA transcription
- Reverse complement generation
- Protein translation
- Open Reading Frame (ORF) detection
- Restriction enzyme analysis (EcoRI, BamHI, HindIII)
- Generate analysis report (.txt)
- Compare GC content between multiple species
- Visualize GC content using Matplotlib

---

## Technologies Used

- Python 3.14
- BioPython
- Matplotlib
- Git
- GitHub

---

## Project Structure

```
DNA-Sequence-Analysis
│
├── data/
│   └── sample.fasta
│
├── images/
│   └── gc_content_comparison.png
│
├── results/
│   └── analysis_report.txt
│
├── src/
│   ├── main.py
│   └── multiple_sequences.py
│
├── README.md
├── LICENSE
└── .gitignore
```

---

## How to Run

Clone the repository:

```bash
git clone https://github.com/bindukanumuri/DNA-Sequence-Analysis.git
```

Install dependencies:

```bash
pip install biopython matplotlib
```

Run the program:

```bash
python src/main.py
```

For multiple FASTA sequences:

```bash
python src/multiple_sequences.py
```

---

## Example Output

The project generates:

- DNA sequence statistics
- GC & AT content
- RNA sequence
- Reverse complement
- Protein sequence
- ORF information
- Restriction enzyme sites
- Analysis report (.txt)
- GC Content Comparison graph

---

## Author

**Bindu Kanumuri**

B.Tech Bioinformatics Student

Interested in Bioinformatics, Machine Learning, Python, and Artificial Intelligence.
