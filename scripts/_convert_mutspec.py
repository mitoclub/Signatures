"""
format mutspec table to such style for SigProfilerExtractor:

Mutation Types      Sample1 Sample2
A[C>A]A             58      74
A[C>A]C             36      66
A[C>A]G             13      12
A[C>A]T             37      64

C>A, C>G, C>T, T>A, T>C, T>G
"""
import os
import random
from functools import reduce, partial

import pandas as pd
import click

PATH_TO_MUTSPEC_192 = "./data/raw/sars-cov-2_ff_mutspec_192comp.csv"
PATH_TO_OUT = "./data/mutspec_96.tsv"
MAX_CHAR_NUM = 15

COLS = ["NucSubst", "RawMutSpec"]
translator = str.maketrans("ATGC", "TACG")


def reformat_mut_style(nuc_subst: str):
    """ replace AAA>ACA style to A[A>C]A style

    RNA mutspec is allowable

    nuc_subst: AAA>ACA format
    """
    assert isinstance(nuc_subst, str) and len(nuc_subst) == 7, "wrong mutation format"
    nuc_subst = nuc_subst.replace("U", "T")
    mutated_nuc = nuc_subst[1]
    if mutated_nuc not in {"C", "T"}:  # see image https://cancer.sanger.ac.uk/signatures/sbs/sbs1/
        nuc_subst = nuc_subst.translate(translator)  # TODO wrong transformation!!!!!!!!!!!!!!, because of reversing absence

    down_nuc = nuc_subst[0]
    up_nuc = nuc_subst[-1]
    mut = nuc_subst[1] + ">" + nuc_subst[-2]
    new_style_subst = down_nuc + "[" + mut + "]" + up_nuc
    return new_style_subst


def process_one_mutspec(path: str, label=None, rounding=True):
    # label = label or "Sample_{:02}".format(random.randint(1, 99))
    label = label or os.path.basename(path)[:MAX_CHAR_NUM].rstrip(".csv").replace(".", "")

    mutspec192 = pd.read_csv(path, usecols=COLS)
    assert mutspec192.shape[0] == 192, "Not all mutations in mutspec!!!"
    mutspec192["MutType"] = mutspec192.NucSubst.apply(reformat_mut_style)

    mutspec96 = mutspec192.groupby("MutType").RawMutSpec.sum().reset_index()
    mutspec96.columns = ["Mutation Types", label]
    if rounding:
        mutspec96[label] = mutspec96[label].round(0).astype(int)
    return mutspec96


def write_mutspec(data: pd.DataFrame, path):
    data.to_csv(path, index=None, sep="\t")


@click.command("converter", help="Convert 192-component mutational spectra to 96-component format")
@click.argument("mutspec_path", nargs=-1, required=True, type=click.Path(True))
@click.option("--round/--no-round", default=False, show_default=True, help="maximum number of signatures to release")
@click.option("--out",  required=True, type=click.Path(), help="path to output mutspec96-sample table")
def main(mutspec_path, round, out):
    mutspec = []
    for path in mutspec_path:
        m = process_one_mutspec(path, rounding=round)
        mutspec.append(m)
    
    mutspec = reduce(partial(pd.merge, on="Mutation Types"), mutspec)
    write_mutspec(mutspec, out)


if __name__ == "__main__":
    main()
