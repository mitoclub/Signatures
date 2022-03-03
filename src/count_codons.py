import json
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

PATH_TO_GENOME = "./data/external/human_genome/GCF_000001405.25_GRCh37.p13_genomic.fna"
PATH_TO_JSON_OUT = "./data/interim/codon_counts_GRCh37.json"

NUCL_SET = set("ACGTacgt")


def is_appropriate_codon(codon: str) -> bool:
    return len(set(codon).difference(NUCL_SET)) == 0


def main():
    codon_counts = defaultdict(int)
    fasta = SeqIO.parse(PATH_TO_GENOME, "fasta")
    rec: SeqRecord = None
    for rec in fasta:
        print("Processing...", rec.description, file=sys.stderr)
        seq = str(rec.seq)
        # iterate over codons with window=1
        for i in range(len(seq) - 2):
            codon = seq[i: i + 3]
            if is_appropriate_codon(codon):
                codon_counts[codon] += 1

    print(codon_counts, file=sys.stderr)
    with open(PATH_TO_JSON_OUT, "w") as fout:
        json.dump(codon_counts, fout)


if __name__ == "__main__":
    main()
