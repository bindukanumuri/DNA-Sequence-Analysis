from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch

print("DNA Sequence Analysis")
print("=" * 40)

sequence = input("Enter a DNA sequence: ").upper()

dna = Seq(sequence)

valid_bases = {"A", "T", "G", "C"}

for base in dna:
    if base not in valid_bases:
        print(f"Invalid DNA sequence! '{base}' is not a valid nucleotide.")
        exit()

print("\nDNA Sequence:", dna)
print("Length:", len(dna))

print("\nNucleotide Count")
print("----------------")

print("A:", dna.count("A"))
print("T:", dna.count("T"))
print("G:", dna.count("G"))
print("C:", dna.count("C"))

gc_content = ((dna.count("G") + dna.count("C")) / len(dna)) * 100
at_content = ((dna.count("A") + dna.count("T")) / len(dna)) * 100

print(f"\nGC Content : {gc_content:.2f}%")
print(f"AT Content : {at_content:.2f}%")

rna = dna.transcribe()

print("\nRNA Sequence")
print(rna)

reverse_comp = dna.reverse_complement()

print("\nReverse Complement")
print(reverse_comp)

if len(dna) % 3 == 0:
    protein = dna.translate()
    print("\nProtein Sequence")
    print(protein)
else:
    print("\nProtein Sequence")
    print("Cannot translate because sequence length is not a multiple of 3.")

with open("results/analysis_report.txt", "w") as file:

    file.write("DNA Sequence Analysis Report\n")
    file.write("=" * 35 + "\n\n")

    file.write(f"DNA Sequence:\n{dna}\n\n")
    file.write(f"Length: {len(dna)}\n\n")

    file.write("Nucleotide Count\n")
    file.write("----------------\n")
    file.write(f"A: {dna.count('A')}\n")
    file.write(f"T: {dna.count('T')}\n")
    file.write(f"G: {dna.count('G')}\n")
    file.write(f"C: {dna.count('C')}\n\n")

    file.write(f"GC Content: {gc_content:.2f}%\n")
    file.write(f"AT Content: {at_content:.2f}%\n\n")

    file.write(f"RNA Sequence:\n{rna}\n\n")
    file.write(f"Reverse Complement:\n{reverse_comp}\n\n")

    if len(dna) % 3 == 0:
        file.write(f"Protein Sequence:\n{protein}\n")
    else:
        file.write("Protein Sequence:\nCannot translate because sequence length is not a multiple of 3.\n")

print("\nReport saved to results/analysis_report.txt")
# ======================================
# Open Reading Frame (ORF)
# ======================================

print("\nOpen Reading Frames")
print("-------------------")

stop_codons = ["TAA", "TAG", "TGA"]

for i in range(len(dna) - 2):

    codon = dna[i:i+3]

    if codon == "ATG":

        for j in range(i + 3, len(dna) - 2, 3):

            stop = dna[j:j+3]

            if stop in stop_codons:

                orf = dna[i:j+3]

                print(f"Start : {i+1}")
                print(f"Stop  : {j+3}")
                print(f"ORF   : {orf}")
                print("-" * 40)

                break

# ======================================
# Restriction Enzyme Analysis
# ======================================

print("\nRestriction Enzyme Analysis")
print("---------------------------")

enzymes = RestrictionBatch(["EcoRI", "BamHI", "HindIII"])

analysis = enzymes.search(dna)

for enzyme, sites in analysis.items():
    print(f"{enzyme}: {sites}")