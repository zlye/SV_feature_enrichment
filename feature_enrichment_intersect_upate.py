#!/usr/bin/env python
# coding: utf-8

import pybedtools
from pybedtools import BedTool
import pandas as pd
import itertools
import datetime

# Use simulations with pybedtools to look for enrichment of
# certain features

gff_path = "/Users/zoelye/Documents/CNV_215_indicaref/R498_refgenome/R498_IGDBv3_coreset.gff"
data_path = "/Users/zoelye/Documents/CNV_215_indicaref/merged_results/merged_duphold_overlap_filter_tableformat.txt"
fai = "/Users/zoelye/Documents/CNV_215_indicaref/R498_refgenome/R498_Chr.fasta.fai"
outfile = "/Users/zoelye/Documents/CNV_215_indicaref/feature_enrichment/SV_iter5000"

minlen = 50
maxlen = 100000
iter_num = 5000

run_big_seq_set = False
run_all_sv = False
run_svtypes = True


def main():
    sv_dat = filter_raw(data_path, minlen, maxlen)
    print("read in dat")
    chromsizes = parse_fai(fai)
    print("generated fai")

    feat_list, bed_list = make_feat_dict(gff_path, run_big_seq_set, chromsizes)
    bed_sv, result_keys = make_sv_bed(chromsizes, sv_dat)
    if run_all_sv is True:
        run_iter_all_sv(feat_list, bed_list, bed_sv, iter_num, result_keys)
    if run_svtypes is True:
        run_iter_sv_type(feat_list, chromsizes, bed_list,
                         sv_dat, iter_num, result_keys)


def filter_raw(cnv_meta, minlen, maxlen):
    # read in data subset
    cnv_data = pd.read_table(data_path)
    cnv_meta = cnv_data.iloc[:, 0:5]
    # FILTERING FOR SIZE: DEL DUP INV, REMOVE INS
    cnv_meta["LEN"] = cnv_meta["END"] - cnv_meta["POS"]

    tmp = cnv_meta[~((cnv_meta["LEN"] < minlen) & (
        cnv_meta["SVTYPE"] == "DEL")) & (cnv_meta["SVTYPE"] != "INS")]

    tmp2 = tmp[~((tmp["LEN"] < minlen) & (tmp["SVTYPE"] == "DUP")) &
               ~((tmp["LEN"] < minlen) & (tmp["SVTYPE"] == "INV"))]
    cnv_meta_size = tmp2[tmp2['LEN'] < maxlen].copy()
    print("svtype counts: ", tmp2.groupby('SVTYPE').count())
    return(cnv_meta_size)


def make_feat_dict(gff_path, run_big_seq_set, chromsizes):
    # read in bed files
    # make dict of feature specific bedfiles
    print("making feat dict")
    gff = pybedtools.BedTool(gff_path).remove_invalid().saveas()
    gff_gene = gff.filter(lambda x: x[2] == 'gene').saveas()
    gff_exon = gff.filter(lambda x: x[2] == 'exon').saveas()
    coding = gff.filter(lambda x: x[2] == 'CDS').saveas()
    gff_three_utr = gff.filter(lambda x: x[2] == 'three_prime_utr').saveas()
    gff_five_utr = gff.filter(lambda x: x[2] == 'five_prime_utr').saveas()
    gff_chr = gff.filter(lambda x: x[2] == 'chromosome').saveas()

    # make the introns - ranges intersect a gene but not coding or UTR
    intwindows = gff_gene.subtract(gff_exon)
    introns = intwindows.intersect(gff_three_utr, v=True)\
                        .intersect(gff_five_utr, v=True).saveas()
    # make intergenic - chr windows minus gene windows
    intergenic = gff_chr.subtract(gff_gene).saveas()
    # make upstream
    upstream = gff_gene.set_chromsizes(chromsizes).flank(l=2000, r=0, s=True)

    # make dicts
    if run_big_seq_set is True:
        print("making big seq list")
        feat_list = ['gene', 'exon', 'intron',
                     'three_prime_utr', 'five_prime_utr', 'intergenic', 'coding', 'upstream']
        bed_list = [gff_gene, gff_exon, introns,
                    gff_three_utr, gff_five_utr, intergenic, coding, upstream]

    else:
        print("making minimal feat list")
        feat_list = ['gene', 'intron',
                     'intergenic', 'coding', 'upstream']
        bed_list = [gff_gene, introns,
                    intergenic, coding, upstream]
    return feat_list, bed_list


def parse_fai(fai):
    """
    parse the fai file in to a dict
    """
    chrom_dict = {}
    with open(fai) as fai_file:
        for line in fai_file:
            line_elems = line.strip().split()
            chrom = line_elems[0]
            chrom_length = int(line_elems[1])
            chrom_dict[chrom] = (0, chrom_length)
    return(chrom_dict)


def make_sv_bed(chromsizes, sv_dat):
    # make bed_sv for all
    bed_sv = pybedtools.BedTool.from_dataframe(
        sv_dat[['CHROM', 'POS', 'END', 'SVTYPE', 'ID']]).set_chromsizes(chromsizes)
    result_keys = ['self', 'other', 'actual', 'frac randomized above actual',
                   'frac randomized below actual', 'median randomized',
                   'normalized', 'percentile']
    return(bed_sv, result_keys)


def run_iter_all_sv(feat_list, bed_list, bed_sv, iter_num, result_keys, outfile):
    # Run overlap iterations FOR ALL SVs
    all_result = {}
    for feature, bed in zip(feat_list, bed_list):
        # run iteraction and write output
        result_dict = bed_sv.randomstats(bed,
                                         iterations=iter_num,
                                         shuffle_kwargs={'chrom': True},
                                         debug=True)
        result_list = [result_dict[x] for x in result_keys]
        all_result[feature] = result_list

        result_df = pd.DataFrame.from_dict(all_result, orient='index',
                                           columns=result_keys)
        result_df['feature'] = feat_list
        # enrichment score: median of the randomized overlap, divided bu actual
        result_df['enrichment_score'] = result_df['actual'] / \
            result_df['median randomized']
        # write output
        out = "_".join([outfile, "sv_enrich_result.csv"])
        result_df.to_csv(out, header=True, index=False)


def run_iter_sv_type(feat_list, chromsizes, bed_list, sv_dat, iter_num, result_keys):
    # MAKE LIST OF PAIRS - SV types and features
    # ITERATIONS FOR PAIRS = SEQ TYPE and SV
    # list of bed files for each sv type
    print('running iterations on all sv types')
    svtype_list = ['BND', 'DEL', 'DUP', 'INV']
    sv_dict = dict(tuple(sv_dat.groupby('SVTYPE')))
    feat_dict = dict(zip(feat_list, bed_list))

    sv_bed_dict = {}
    for i in svtype_list:
        sv_bed_dict[i] = pybedtools.BedTool.from_dataframe(
            sv_dict[i][['CHROM', 'POS', 'END', 'SVTYPE', 'ID']]
        ).set_chromsizes(chromsizes)
    print("bed filter for each SV type generated")

    svtype_feat_pairs = list(itertools.product(svtype_list, feat_list))

    result_dat = {}
    for pair in svtype_feat_pairs:
        sv_bed = sv_bed_dict[pair[0]]
        feat_bed = feat_dict[pair[1]]
        result_dict = sv_bed.randomstats(feat_bed,
                                         iterations=iter_num,
                                         shuffle_kwargs={'chrom': True},
                                         debug=True)
        print("ran random stats for {}".format(pair))
        print(datetime.datetime.now())
        result_list = [result_dict[x] for x in result_keys]
        result_dat[pair] = result_list

    result_df = pd.DataFrame.from_dict(result_dat, orient='index',
                                       columns=result_keys)
    # make data frame of results
    idx_unzip = list(zip(*result_df.index.tolist()))
    result_df['svtype'] = idx_unzip[0]
    result_df['feature'] = idx_unzip[1]
    # enrichment score: median of the randomized overlap, divided bu actual
    result_df['enrichment_score'] = result_df['actual'] / \
        result_df['median randomized']
    print("enrichment scrote calculation")
    # write output
    out = "_".join([outfile, "svtype_enrich_result.csv"])
    result_df.to_csv(out, header=True, index=False)


if __name__ == '__main__':
    main()
