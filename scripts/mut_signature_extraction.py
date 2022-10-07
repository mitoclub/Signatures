import numpy as np
import pandas as pd
from SigProfilerExtractor import sigpro as sig
import click


def unique_py(seqlist):
    seen = set()
    seen_add = seen.add
    return [x for x in seqlist if not (x in seen or seen_add(x))]


def calcIntermutDist2(subs_type, first_chrom_na=False):
    subs_type_processed = subs_type.copy()
    chr_list = unique_py(subs_type['chr'])
    pos_array_im = subs_type['position'].values
    index_orig_df = np.arange(len(subs_type_processed))
    #args_pos_list = np.argsort(pos_array_im)
    args_pos_list = []
    distPrev_list = []
    prevPos_list = []

    for c in chr_list:
        inds_chr = np.where(subs_type['chr'] == c)
        pos_array_im_c = np.sort(pos_array_im[inds_chr])
        index_orig_df[inds_chr] = index_orig_df[inds_chr][np.argsort(
            pos_array_im[inds_chr])]

        if first_chrom_na:
            prevPos_arr_c = np.hstack((np.NAN, pos_array_im_c.flatten()[:-1]))
        else:
            prevPos_arr_c = np.hstack((0, pos_array_im_c.flatten()[:-1]))
        distPrev_arr_c = pos_array_im_c - prevPos_arr_c
        distPrev_arr_c[distPrev_arr_c == 0] = 1
        distPrev_list = np.append(
            distPrev_list, distPrev_arr_c.astype(int)).flatten()
        prevPos_list = np.append(
            prevPos_list, prevPos_arr_c.astype(int)).flatten()
        prevPos_arr_c = []
        distPrev_arr_c = []
    subs_type_processed = subs_type_processed.reindex(
        index_orig_df).reset_index(drop=True)
    subs_type_processed['prevPos'] = prevPos_list
    subs_type_processed['distPrev'] = distPrev_list
    return subs_type_processed


def calcIntermutDist(subs_type, first_chrom_na=False):
    subs_type_processed = pd.DataFrame()
    for c in unique_py(subs_type['chr']):
        subs_type_chrom = subs_type[subs_type['chr']
                                    == c].sort_values('position')
        if first_chrom_na:
            subs_type_chrom['prevPos'] = np.hstack(
                (np.NAN, subs_type_chrom['position'].values.flatten()[:-1]))
        else:
            subs_type_chrom['prevPos'] = np.hstack(
                (0, subs_type_chrom['position'].values.flatten()[:-1]))
        subs_type_chrom['distPrev'] = subs_type_chrom['position'].values - \
            subs_type_chrom['prevPos'].values
        subs_type_processed = subs_type_processed.append(subs_type_chrom)
        subs_type_processed['distPrev'][subs_type_processed['distPrev'] == 0] = 1
    return subs_type_processed


def computeIMD3(chrom_df, chromosome):

    # keep track of partners

    d1 = dict(zip(list(chrom_df['start1']), list(chrom_df['start2'])))
    d2 = dict(zip(list(chrom_df['start2']), list(chrom_df['start1'])))
    d = {**d1, **d2}  # combine dictionaries

    lb = chrom_df.iloc[:, 0:2]  # get chrom1 and start1
    rb = chrom_df.iloc[:, 3:5]  # get chrom2 and start2
    rest = chrom_df.iloc[:, 6:]

    lb = pd.DataFrame(np.concatenate((lb.values, rest.values), axis=1))
    rb = pd.DataFrame(np.concatenate((rb.values, rest.values), axis=1))

    # BREAKPOINTS ARE CONSIDERED INDIVIDUALLY

    #['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'sample', 'svclass', 'size_bin', 'length']
    lb.columns = ['chrom1', 'start1', 'sample',
                  'svclass', 'size_bin', "length"]
    rb.columns = ['chrom2', 'start2', 'sample',
                  'svclass', 'size_bin', "length"]

    chr_lb = lb[lb.chrom1 == chromosome]
    chr_rb = rb[rb.chrom2 == chromosome]
    # print(chr_lb)
    # print(chr_rb)
    chrom_df = pd.DataFrame(np.concatenate(
        (chr_lb.values, chr_rb.values), axis=0))
    chrom_df.columns = ['chrom', 'start',
                        'sample', 'svclass', 'size_bin', "length"]

    # print(chrom_df['chrom'].unique())
    #assert(chrom_df['chrom'].nunique() == 1)

    # sort on last column which is start coordinate
    chrom_df = chrom_df.sort_values(chrom_df.columns[1])  # CHROM, START

    # take care of mirrored translocations
    to_drop = []
    starts = list(chrom_df["start"])
    svtypes = list(chrom_df["svclass"])
    for i, (s, svtype) in enumerate(zip(starts, svtypes)):
        if i+1 < len(starts) and abs(starts[i+1] - s) <= 100 and svtype == "translocation":
            to_drop.append(i)

    chrom_df = chrom_df.drop(to_drop)
    chrom_df = chrom_df.sort_values(chrom_df.columns[1])

    coords = list(chrom_df[chrom_df.columns[1]])
    svtype = list(chrom_df.svclass)

    chrom_inter_distances = []

    # defined as the number of base pairs from one rearrangement breakpoint to the one closest to it that is not it's partner
    for i in range(1, len(coords)-1):

        j = i-1
        k = i+1
        # check if previous breakpoint is partner of this breakpoint, if it is, avoid it
        while j >= 0 and coords[j] == d[coords[i]]:
            j = j-1
        while k < len(coords) and coords[k] == d[coords[i]]:
            k = k+1
        if j >= 0 and k < len(coords):
            if coords[i] - coords[j] == 0:
                dist = coords[k] - coords[i]
            elif coords[k] - coords[i] == 0:
                dist = coords[i] - coords[j]
            else:
                dist = min(coords[i] - coords[j], coords[k] - coords[i])
        elif j < 0:
            dist = coords[k] - coords[i]
        else:
            dist = coords[i] - coords[j]

        if dist == 0 and svtype[i] == "translocation":
            print(coords[j], coords[i], coords[k], dist)
            # print(len(coords))
        chrom_inter_distances.append(dist)
        if dist == 1:
            print(coords[j], coords[i], coords[k], svtype[i])

    # now we take care of the edge cases of the first and last breakpoint

    if coords[1] == d[coords[0]]:
        first_dist = coords[2] - coords[0]
    else:
        first_dist = coords[1] - coords[0]

    if coords[-2] == d[coords[-1]]:
        last_dist = coords[-1] - coords[-3]
    else:
        last_dist = coords[-1] - coords[-2]

    chrom_inter_distances = [first_dist] + chrom_inter_distances
    chrom_inter_distances.append(last_dist)
    chrom_df['IMD'] = chrom_inter_distances

    # INTERLEAVED VS NESTED CONFIGURATION
    configuration = ['interleaved' for i in range(len(coords))]
    for i in range(1, len(coords)):
        j = i-1
        # check if previous breakpoint is partner of this breakpoint, if it is, avoid it
        while coords[j] == d[coords[i]] and not (d[coords[i]] < max(d[coords[j]], coords[j]) and coords[i] < max(d[coords[j]], coords[j]) and d[coords[i]] > min(d[coords[j]], coords[j]) and coords[i] > min(d[coords[j]], coords[j])):
            j = j-1
        if j >= 0:  # determine if we have a nested or interleaved configuration
            if d[coords[i]] < max(d[coords[j]], coords[j]) and coords[i] < max(d[coords[j]], coords[j]) and d[coords[i]] > min(d[coords[j]], coords[j]) and coords[i] > min(d[coords[j]], coords[j]):
                configuration[i] = "nested"

    chrom_df["Configuration"] = configuration
    return chrom_df


@click.command("extractor", help="Pass path to mutspec table in SigProfiler format and path to output directory. Script release signatures from mutspec and decompose it with COSMIC database")
@click.argument("input_mutspec", required=True, type=click.Path(exists=True))
@click.argument("outdir", required=True, type=click.Path(exists=False))
@click.option("-m", "--max_signatures", default=5, show_default=True, type=int, help="maximum number of signatures to release")
@click.option("-t", "--threads", default=-1, show_default=True, type=int, help="number of threads to use")
def main(input_mutspec, outdir, max_signatures, threads):
    sig.sigProfilerExtractor(
        "matrix",
        outdir, input_mutspec,
        minimum_signatures=1,
        maximum_signatures=max_signatures,
        cpu=threads,
        gpu=False,
    )


if __name__ == "__main__":
    main()
