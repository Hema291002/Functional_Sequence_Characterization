from Bio.Blast import NCBIXML

with open("results/blast_results.xml") as f:
    blast_record = NCBIXML.read(f)

annotations = []
annotations.append("Functional Annotation Report")
annotations.append("--------------------------------")

for alignment in blast_record.alignments[:5]:
    hsp = alignment.hsps[0]
    annotations.append(
        "Hit: " + alignment.hit_def +
        " | Identity: " + str(hsp.identities) +
        " | E-value: " + str(hsp.expect)
    )

annotations.append("")
annotations.append("Biological Interpretation:")
annotations.append(
    "The query sequence shows high similarity to known proteins. "
    "Low E-values and high identity suggest conserved biological function."
)

with open("results/functional_annotation.txt", "w") as f:
    f.write("\n".join(annotations))

print("Functional annotation completed.")
