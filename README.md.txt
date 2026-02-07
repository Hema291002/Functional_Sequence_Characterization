# Functional Sequence Characterization Pipeline

## Project Title
In-silico Identification and Functional Characterization of a Target Gene/Protein  
Using Sequence Analysis and Homology-Based Annotation

---

## Core Research Question
**Given an unknown or hypothetical biological sequence, how can computational analysis be used to assess its quality, identify homologous sequences, and predict its biological function?**

This question drives the entire pipeline implemented in this project.

---

## Project Overview
This project implements a **single continuous bioinformatics pipeline** that starts with an unknown biological sequence and progressively analyzes it using sequence-based and homology-based approaches to infer its biological function.

The pipeline is implemented using **Biopython** and follows real bioinformatics research logic.

---

## Pipeline Workflow

Sequence Retrieval
↓
Sequence Quality Analysis
↓
Sequence Filtering & Validation
↓
Homology Search (BLAST)
↓
Homology-Based Functional Annotation
↓
Biological Interpretation


Each step depends on the output of the previous step.

---

## Step-wise Description of the Pipeline

Before running the pipeline, create an empty `results/` directory:

### Step 1: Sequence Retrieval
An unknown or hypothetical biological sequence is provided in FASTA format.

**Input file:**  
`data/input_sequence.fasta`

The sequence may represent a hypothetical protein or an uncharacterized gene product.

---

### Step 2: Sequence Quality Analysis
The sequence is first evaluated to determine whether it is biologically reasonable.

Implemented in:  
`analysis/sequence_qc.py`

**Analyses performed:**
- Sequence length calculation
- Identification of sequence type (DNA or protein)
- GC content calculation for DNA sequences
- Amino acid composition analysis for protein sequences

---

### Step 3: Sequence Filtering & Validation
Based on biological criteria, a decision is made on whether the sequence is suitable for downstream analysis.

**Filtering rule applied:**
- Sequences shorter than 50 residues are rejected
- Longer sequences are accepted for homology analysis

**Output file:**  
`results/qc_summary.txt`

This step simulates real research decision-making prior to functional analysis.

---

### Step 4: Homology Search (BLAST)
Homology search is performed using **NCBI BLAST through Biopython**.

Implemented in:  
`analysis/run_blast.py`

**Method:**
- BLASTP is used for protein-protein similarity search
- The query sequence is compared against the NCBI non-redundant (nr) protein database
- The BLAST search is submitted online using Biopython’s `NCBIWWW` interface

**Output file:**  
`results/blast_results.xml`

Homology search is the central step of functional prediction.

---

### Step 5: Homology-Based Functional Annotation
Functional annotation is inferred from the BLAST results.

Implemented in:  
`analysis/homology_analysis.py`

**Approach:**
- BLAST XML output is parsed using `NCBIXML`
- Top homologous hits are extracted
- Sequence identity and E-value are used as evidence for functional similarity

**Output file:**  
`results/functional_annotation.txt`

---

### Step 6: Biological Interpretation
The final step involves interpreting the computational evidence to infer biological function.

Low E-values and high sequence identity indicate conserved function, allowing prediction of the likely role of the query sequence.

This step transforms raw computational output into a research conclusion.

---

## Project Structure

Functional_Sequence_Characterization/
├── analysis/
│ ├── sequence_qc.py
│ ├── run_blast.py
│ └── homology_analysis.py
├── data/
│ └── input_sequence.fasta
├── results/
│ ├── qc_summary.txt
│ ├── blast_results.xml
│ ├── functional_annotation.txt
│ └── query.fasta
├── .gitignore
└── README.md


---

## Tools and Libraries Used
- Python 3
- Biopython
  - Seq
  - SeqIO
  - SeqRecord
  - NCBIWWW
  - NCBIXML

---

## How to Run the Pipeline

From the project root directory:

```bash
python analysis/sequence_qc.py
python analysis/run_blast.py
python analysis/homology_analysis.py
Note: An active internet connection is required for the BLAST step.

Final Outcome
By completing this pipeline, the user can state:

“An unknown biological sequence was analyzed for quality, filtered based on biological criteria, compared against known sequences using BLAST, and functionally annotated based on homology evidence.”

This project reflects real bioinformatics research methodology.

Author
Hema 

