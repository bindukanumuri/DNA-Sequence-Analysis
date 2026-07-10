from Bio import SeqIO
import matplotlib.pyplot as plt

report = open("results/analysis_report.txt", "w")
sequence_names = []
gc_values = [] 

print("DNA Sequence Analysis")
print("=" * 40)

report.write("DNA Sequence Analysis\n")
report.write("=" * 40 + "\n\n")

for record in SeqIO.parse("data/sample.fasta", "fasta"):

    dna = record.seq

    print(f"\nSequence ID : {record.id}")
    print(f"DNA         : {dna}")
    print(f"Length      : {len(dna)}")

    report.write(f"Sequence ID : {record.id}\n")
    report.write(f"DNA         : {dna}\n")
    report.write(f"Length      : {len(dna)}\n")

    print("\nNucleotide Count")
    print("----------------")

    report.write("\nNucleotide Count\n")
    report.write("----------------\n")

    print("A:", dna.count("A"))
    print("T:", dna.count("T"))
    print("G:", dna.count("G"))
    print("C:", dna.count("C"))

    report.write(f"A: {dna.count('A')}\n")
    report.write(f"T: {dna.count('T')}\n")
    report.write(f"G: {dna.count('G')}\n")
    report.write(f"C: {dna.count('C')}\n")

    gc = ((dna.count("G") + dna.count("C")) / len(dna)) * 100
    at = ((dna.count("A") + dna.count("T")) / len(dna)) * 100

    sequence_names.append(record.id)
    gc_values.append(gc)

    print(f"\nGC Content : {gc:.2f}%")
    print(f"AT Content : {at:.2f}%")

    report.write(f"\nGC Content : {gc:.2f}%\n")
    report.write(f"AT Content : {at:.2f}%\n")

    rna = dna.transcribe()
    reverse = dna.reverse_complement()

    report.write(f"\nRNA Sequence:\n{rna}\n")
    report.write(f"\nReverse Complement:\n{reverse}\n")

    if len(dna) % 3 == 0:
        protein = dna.translate()
        report.write(f"\nProtein Sequence:\n{protein}\n")
    else:
        report.write("\nProtein Sequence:\nCannot translate (length not multiple of 3).\n")

    print("-" * 50)
    report.write("\n" + "-" * 50 + "\n\n")

report.close()

print("\nAnalysis report saved to results/analysis_report.txt")

plt.figure(figsize=(6,4))

plt.bar(sequence_names, gc_values)

plt.title("GC Content Comparison")
plt.xlabel("Species")
plt.ylabel("GC Content (%)")

plt.savefig("images/gc_content_comparison.png")

plt.show()