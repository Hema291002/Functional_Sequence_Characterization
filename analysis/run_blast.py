from Bio import SeqIO
from Bio.Blast import NCBIWWW

# Run from project root directory
# Internet connection required

record = SeqIO.read("data/input_sequence.fasta", "fasta")

print("Submitting BLASTP search using Biopython...")

result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=record.seq
)

with open("results/blast_results.xml", "w") as f:
    f.write(result_handle.read())

print("BLAST completed. Results saved to results/blast_results.xml")
