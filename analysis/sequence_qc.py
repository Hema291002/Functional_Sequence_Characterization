from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Run from project root directory

record = SeqIO.read("data/input_sequence.fasta", "fasta")
seq = record.seq
length = len(seq)

dna_bases = set("ATGC")
is_dna = all(base in dna_bases for base in seq)

results = []
results.append("Sequence ID: " + record.id)
results.append("Sequence Length: " + str(length))

if is_dna:
    g = seq.count("G")
    c = seq.count("C")
    gc = ((g + c) / length) * 100
    results.append("Sequence Type: DNA")
    results.append("GC Content: " + str(round(gc, 2)) + "%")
else:
    results.append("Sequence Type: Protein")
    results.append("Amino Acid Composition:")
    for aa in sorted(set(seq)):
        results.append(aa + ": " + str(seq.count(aa)))

# Filtering decision
if length < 50:
    results.append("Filtering Decision: Sequence rejected (too short)")
else:
    results.append("Filtering Decision: Sequence accepted for homology analysis")

with open("results/qc_summary.txt", "w") as f:
    f.write("\n".join(results))

print("Sequence quality analysis and filtering completed.")
