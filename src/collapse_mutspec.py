import json
import os
import random
from collections import defaultdict
from functools import partial, reduce

import click
import numpy as np
import pandas as pd

PATH_TO_HUMAN_COUNTS = "./data/codon_counts_GRCh37.json"
COLS = ["NucSubst", "RawMutSpec"]
translator = str.maketrans("ATGC", "TACG")


def load_triplet_counts(path: str, normalize=False) -> pd.Series:
    """ read and collapse raw trinucleotide counts """
    considerable_nucs = {"C", "T"}
    with open(path) as fin:
        counts = json.load(fin)

    new_counts = defaultdict(int)
    for trinuc, num in counts.items():
        standart_trinuc = trinuc.upper()
        if len(set(standart_trinuc).difference("ATGC")) == 0:
            if standart_trinuc[1] not in considerable_nucs:
                standart_trinuc = standart_trinuc.translate(translator)
            new_counts[standart_trinuc] += num

    new_counts = pd.Series(new_counts)
    if normalize:
        new_counts = new_counts / new_counts.median()
    return new_counts.sort_index()


def reformat_mut_style(nuc_subst: str):
    """ replace AAA>ACA style to A[A>C]A style

    RNA mutspec is allowable

    nuc_subst: AAA>ACA format
    """
    assert isinstance(nuc_subst, str) and len(
        nuc_subst) == 7, "wrong mutation format"
    nuc_subst = nuc_subst.replace("U", "T")
    revcomp = False
    # see image https://cancer.sanger.ac.uk/signatures/sbs/sbs1/
    if nuc_subst[1] not in {"C", "T"}:
        nuc_subst = nuc_subst.translate(translator)
        revcomp = True

    down_nuc, up_nuc = nuc_subst[0], nuc_subst[-1]
    N1, N2 = nuc_subst[1], nuc_subst[-2]

    if revcomp:
        sbs96 = up_nuc + "[" + N1 + ">" + N2 + "]" + down_nuc
    else:
        sbs96 = down_nuc + "[" + N1 + ">" + N2 + "]" + up_nuc
    return sbs96


def collapse_mutspec(mutspec192: pd.DataFrame):
    assert mutspec192.shape[0] == 192, "Not all mutations in mutspec!!!"
    if mutspec192["NucSubst"].str.contains("^\w{3}>\w{3}$").sum() == 192:
        # collapsing: AGG>ATG -> A[G>T]G -> C[C>A]T
        mutspec192["Mutation Types"] = mutspec192.NucSubst.apply(reformat_mut_style)
        mutspec192.drop("NucSubst", axis=1, inplace=True)
        mutspec96 = mutspec192.groupby("Mutation Types").sum().reset_index()
    else:
        raise NotImplementedError("Such MutSpec type is not allovable")

    cols = ["Mutation Types"] + list(mutspec96.columns.drop("Mutation Types"))
    mutspec96 = mutspec96[cols]
    return mutspec96


def renormalize(mutspec, inplace=False, scale=True):
    if not inplace:
        mutspec = mutspec.copy()
    mutspec["Context"] = mutspec["Mutation Types"].str.extract(
        "(\w)\[(\w)>\w\](\w)").apply(lambda x: "".join(x), axis=1)

    human_triplet_counts = load_triplet_counts(PATH_TO_HUMAN_COUNTS, True)
    triplet_freqs = human_triplet_counts.reset_index(name="Freq")\
        .rename({"index": "Context"}, axis=1)
    ext_mutspec = pd.merge(mutspec, triplet_freqs, on="Context")\
        .sort_values("Mutation Types")

    additional_cols = {"Context", "Freq", "Mutation Types"}
    for c in ext_mutspec.columns:
        if c in additional_cols or ext_mutspec[c].max() <= 0:
            continue
        new_mutspec = ext_mutspec[c] * ext_mutspec["Freq"]
        if scale:
            mina = new_mutspec[new_mutspec > 0].min()
            new_mutspec = np.round(new_mutspec / mina).astype(int)
        ext_mutspec[c] = new_mutspec
    ext_mutspec.drop(["Context", "Freq"], axis=1, inplace=True)
    return ext_mutspec


def write_mutspec(data: pd.DataFrame, path):
    data.to_csv(path, index=None, sep="\t")


@click.command("converter", help="Convert 192-component mutational spectra to 96-component format")
@click.argument("mutspec_path", required=True, type=click.Path(True))
@click.option("--out",  required=True, type=click.Path(), help="path to output mutspec96-sample table (tsv)")
@click.option("--scale/--no-scale", default=True, show_default=True, help="Divide by minimal value in sample mutspec or not")
def main(mutspec_path, out, scale):
    mutspec192 = pd.read_csv(mutspec_path)
    mutspec96 = collapse_mutspec(mutspec192)
    mutspec96human = renormalize(mutspec96, scale=scale)
    write_mutspec(mutspec96human, out)


if __name__ == "__main__":
    # main("./data/raw/birds192.csv", "/tmp/file.tsv", True)
    main()
