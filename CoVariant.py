#!/bin/python3
import pandas as pd
import os
import fnmatch
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Sample_List", help="A tab delineated file with each line containing sample names. Last line is empty.")
parser.add_argument("Virus", help="MHV (AY910861.1), MERS (JX869059.2), or SARS2 (MT020881.1)")
parser.add_argument("Working_Directory", help="Path to directory containing data to align.")
parser.add_argument("Experiment", help="Experiment name.")
parser.add_argument("--freq", help="Variant frequency cutoff for filtering. Decimal between 0 and 1.")
parser.add_argument("--file_tag", help="File naming tag can denote filters used for variant isolation.")
args = parser.parse_args()

sample_list = [line.rstrip('\n') for line in open(str(args.Sample_List))]
# For virus, choose "MHV", "MERS", or "SARS2"
virus = str(args.Virus)
wd = str(args.Working_Directory)
od = wd
exp = str(args.Experiment)
if args.file_tag:
    tag = "_" + str(args.file_tag)
else:
    tag = ""
if args.freq:
    freq_cutoff = float(args.freq)
else:
    freq_cutoff = 0
report = pd.DataFrame(columns=['sample',
                               'unique_variants',
                               'variant_nts',
                               "total_nts",
                               "transition_nts",
                               "transversion_nts",
                               "AtoG_nts",
                               "GtoA_nts",
                               "CtoT_nts",
                               "TtoC_nts",
                               "AtoT_nts",
                               "TtoA_nts",
                               "AtoC_nts",
                               "CtoA_nts",
                               "CtoG_nts",
                               "GtoC_nts",
                               "GtoT_nts",
                               "TtoG_nts",
                               "mutation_freq",
                               "transition_freq",
                               "transversion_freq",
                               "AtoG_freq",
                               "GtoA_freq",
                               "CtoT_freq",
                               "TtoC_freq",
                               "AtoT_freq",
                               "TtoA_freq",
                               "AtoC_freq",
                               "CtoA_freq",
                               "CtoG_freq",
                               "GtoC_freq",
                               "GtoT_freq",
                               "TtoG_freq"
                               ])
report['sample'] = sample_list

for file in os.listdir(wd):
    if fnmatch.fnmatch(file, "*_coverage.txt"):
        sample_name = file.split("_")[0]
        depth = pd.read_csv(wd + file, sep="\t", header = 0)
        total_depth = sum(depth['Coverage'])
        report.loc[report['sample'] == sample_name, ['total_nts']] = total_depth

for file in os.listdir(wd):
    if fnmatch.fnmatch(file, "*.vcf"):
        sample_name = str(file.split(".")[0])
        vcf = pd.read_csv(wd + file, skiprows=18, sep="\t", header=0, index_col=False, names=["genome",
                                                                                              "position",
                                                                                              "ID",
                                                                                              "reference",
                                                                                              "variant",
                                                                                              "qual",
                                                                                              "filter",
                                                                                              "info"])
        vcf[['raw_depth',
             'frequency',
             'strand_bias',
             'DP4']] = vcf['info'].apply(lambda x: pd.Series(x.split(';')))
        vcf['raw_depth'] = vcf['raw_depth'].str[3:]
        vcf['frequency'] = vcf['frequency'].str[3:]
        vcf['DP4'] = vcf['DP4'].str[4:]
        vcf[['ref_f_count', 'ref_r_count', 'variant_f_count', 'variant_r_count']] = vcf['DP4'].apply(lambda x: pd.Series(x.split(',')))
        vcf = vcf.drop(columns=['ID', "qual", "filter", 'info', 'strand_bias', 'DP4'])
        vcf = vcf[['genome', 'position', 'reference', 'variant', 'frequency', 'raw_depth', 'ref_f_count', 'ref_r_count', 'variant_f_count', 'variant_r_count']]
        vcf['position'] = pd.to_numeric(vcf['position'])
        vcf['frequency'] = pd.to_numeric(vcf['frequency'])
        vcf['raw_depth'] = pd.to_numeric(vcf['raw_depth'])
        vcf['ref_f_count'] = pd.to_numeric(vcf['ref_f_count'])
        vcf['ref_r_count'] = pd.to_numeric(vcf['ref_r_count'])
        vcf['variant_f_count'] = pd.to_numeric(vcf['variant_f_count'])
        vcf['variant_r_count'] = pd.to_numeric(vcf['variant_r_count'])
        vcf['variant_total'] = vcf['variant_f_count'] + vcf['variant_r_count']
        vcf = vcf[vcf['frequency'] >= freq_cutoff]
        def get_variant_type(reference, variant):
            type = ""
            if reference == "A" and variant == "G":
                type = "transition"
            elif reference == "G" and variant == "A":
                type = "transition"
            elif reference == "C" and variant == "T":
                type = "transition"
            elif reference == "T" and variant == "C":
                type = "transition"
            else:
                type = "transversion"
            return type
        def get_SARS2_gene(position):
            gene = ""
            if position > 0 and position < 265:
                gene = "5UTR"
            elif position > 265 and position < 806:
                gene = "nsp1"
            elif position > 805 and position < 2720:
                gene = "nsp2"
            elif position > 2719 and position < 8555:
                gene = "nsp3"
            elif position > 8554 and position < 10055:
                gene = "nsp4"
            elif position > 10054 and position < 10973:
                gene = "nsp5"
            elif position > 10972 and position < 11843:
                gene = "nsp6"
            elif position > 11842 and position < 12092:
                gene = "nsp7"
            elif position > 12091 and position < 12686:
                gene = "nsp8"
            elif position > 12685 and position < 13025:
                gene = "nsp9"
            elif position > 13024 and position < 13442:
                gene = "nsp10"
            elif position > 13441 and position < 16237:
                gene = "nsp12"
            elif position > 16236 and position < 18040:
                gene = "nsp13"
            elif position > 18039 and position < 19621:
                gene = "nsp14"
            elif position > 19620 and position < 20659:
                gene = "nsp15"
            elif position > 20658 and position < 21553:
                gene = "nsp16"
            elif position > 21562 and position < 25385:
                gene = "S protein"
            elif position > 25392 and position < 26221:
                gene = "ORF3a"
            elif position > 26244 and position < 26473:
                gene = "E protein"
            elif position > 26522 and position < 27192:
                gene = "M protein"
            elif position > 27201 and position < 27388:
                gene = "ORF6"
            elif position > 27393 and position < 27888:
                gene = "ORF7ab"
            elif position > 27893 and position < 28260:
                gene = "ORF8"
            elif position > 28273 and position < 29534:
                gene = "N protein"
            elif position > 29557 and position < 29675:
                gene = "ORF10"
            elif position > 29674:
                gene = "3UTR"
            else:
                gene = "unknown"
            return gene
        def get_MHV_gene(position):
            gene = ""
            if position > 0 and position < 210:
                gene = "5UTR"
            elif position > 209 and position < 951:
                gene = "nsp1"
            elif position > 950 and position < 2706:
                gene = "nsp2"
            elif position > 2705 and position < 9633:
                gene = "nsp3"
            elif position > 9632 and position < 10209:
                gene = "nsp4"
            elif position > 10208 and position < 11118:
                gene = "nsp5"
            elif position > 11117 and position < 11979:
                gene = "nsp6"
            elif position > 11978 and position < 12246:
                gene = "nsp7"
            elif position > 12245 and position < 12837:
                gene = "nsp8"
            elif position > 12836 and position < 13167:
                gene = "nsp9"
            elif position > 13166 and position < 13578:
                gene = "nsp10"
            elif position > 13577 and position < 16361:
                gene = "nsp12"
            elif position > 16360 and position < 18161:
                gene = "nsp13"
            elif position > 18160 and position < 19724:
                gene = "nsp14"
            elif position > 19723 and position < 20846:
                gene = "nsp15"
            elif position > 20845 and position < 21743:
                gene = "nsp16"
            elif position > 21744 and position < 21754:
                gene = "TRS-2"
            elif position > 21770 and position < 22557:
                gene = "ORF2a"
            elif position > 22601 and position < 23922:
                gene = "HE"
            elif position > 23919 and position < 23929:
                gene = "TRS-3"
            elif position > 23928 and position < 27904:
                gene = "S protein"
            elif position > 27932 and position < 27942:
                gene = "TRS-4"
            elif position > 27992 and position < 28053:
                gene = "ORF4a"
            elif position > 28057 and position < 28379:
                gene = "ORF4b"
            elif position > 28315 and position < 28325:
                gene = "TRS-5"
            elif position > 28374 and position < 28714:
                gene = "ORF5a"
            elif position > 28705 and position < 28957:
                gene = "E protein"
            elif position > 28955 and position < 28965:
                gene = "TRS-6"
            elif position > 28967 and position < 29655:
                gene = "M protein"
            elif position > 29652 and position < 29662:
                gene = "TRS-7"
            elif position > 29668 and position < 31032:
                gene = "N protein"
            elif position > 31033:
                gene = "3UTR"
            else:
                gene = "unknown"
            return gene
        vcf['variant_type'] = vcf[['reference', 'variant']].apply(lambda x: get_variant_type(*x), axis=1)
        if (virus == "SARS2"):
            print("Using SARS-CoV-2 annotations corresponding to MT020881.1.")
            vcf['gene'] = vcf['position'].apply(lambda x: get_SARS2_gene(x))
        elif (virus == "MHV"):
            print("Using MHV annotations corresponding to AY910861.1.")
            vcf['gene'] = vcf['position'].apply(lambda x: get_MHV_gene(x))
        # elif (virus == "MERS"):
        #     vcf['gene'] = vcf['position'].apply(lambda x: get_MERS_gene(x))
        else:
            print("No virus gene annotations available for that virus! Please check that you have the correct virus specified. Otherwise, contact your developer to input annotations.")
        variant_nts = vcf['variant_total'].sum()
        transition_nts = vcf.loc[vcf['variant_type'] == "transition", 'variant_total'].sum()
        transversion_nts = vcf.loc[vcf['variant_type'] == "transversion", 'variant_total'].sum()
        AtoG_nts = vcf.loc[((vcf['reference'] == "A") & (vcf['variant'] == "G")), 'variant_total'].sum()
        GtoA_nts = vcf.loc[((vcf['reference'] == "G") & (vcf['variant'] == "A")), 'variant_total'].sum()
        CtoT_nts = vcf.loc[((vcf['reference'] == "C") & (vcf['variant'] == "T")), 'variant_total'].sum()
        TtoC_nts = vcf.loc[((vcf['reference'] == "T") & (vcf['variant'] == "C")), 'variant_total'].sum()
        AtoC_nts = vcf.loc[((vcf['reference'] == "A") & (vcf['variant'] == "C")), 'variant_total'].sum()
        CtoA_nts = vcf.loc[((vcf['reference'] == "C") & (vcf['variant'] == "A")), 'variant_total'].sum()
        AtoT_nts = vcf.loc[((vcf['reference'] == "A") & (vcf['variant'] == "T")), 'variant_total'].sum()
        TtoA_nts = vcf.loc[((vcf['reference'] == "T") & (vcf['variant'] == "A")), 'variant_total'].sum()
        CtoG_nts = vcf.loc[((vcf['reference'] == "C") & (vcf['variant'] == "G")), 'variant_total'].sum()
        GtoC_nts = vcf.loc[((vcf['reference'] == "G") & (vcf['variant'] == "C")), 'variant_total'].sum()
        GtoT_nts = vcf.loc[((vcf['reference'] == "G") & (vcf['variant'] == "T")), 'variant_total'].sum()
        TtoG_nts = vcf.loc[((vcf['reference'] == "T") & (vcf['variant'] == "G")), 'variant_total'].sum()
        unique_variants = len(vcf)
        report.loc[report['sample'] == sample_name, ['unique_variants']] = unique_variants
        report.loc[report['sample'] == sample_name, ['variant_nts']] = variant_nts
        report.loc[report['sample'] == sample_name, ['transition_nts']] = transition_nts
        report.loc[report['sample'] == sample_name, ['transversion_nts']] = transversion_nts
        report.loc[report['sample'] == sample_name, ['AtoG_nts']] = AtoG_nts
        report.loc[report['sample'] == sample_name, ['GtoA_nts']] = GtoA_nts
        report.loc[report['sample'] == sample_name, ['CtoT_nts']] = CtoT_nts
        report.loc[report['sample'] == sample_name, ['TtoC_nts']] = TtoC_nts
        report.loc[report['sample'] == sample_name, ['AtoC_nts']] = AtoC_nts
        report.loc[report['sample'] == sample_name, ['CtoA_nts']] = CtoA_nts
        report.loc[report['sample'] == sample_name, ['AtoT_nts']] = AtoT_nts
        report.loc[report['sample'] == sample_name, ['TtoA_nts']] = TtoA_nts
        report.loc[report['sample'] == sample_name, ['CtoG_nts']] = CtoG_nts
        report.loc[report['sample'] == sample_name, ['GtoC_nts']] = GtoC_nts
        report.loc[report['sample'] == sample_name, ['GtoT_nts']] = GtoT_nts
        report.loc[report['sample'] == sample_name, ['TtoG_nts']] = TtoG_nts
        vcf.to_csv(od + sample_name + tag + "_variants.txt", sep="\t", index=False)
report['mutation_freq'] = (report['variant_nts'] / report['total_nts'])
report['transition_freq'] = report['transition_nts'] / report['total_nts']
report['transversion_freq'] = report['transversion_nts'] / report['total_nts']
report['AtoG_freq'] = report['AtoG_nts'] / report['total_nts']
report['GtoA_freq'] = report['GtoA_nts'] / report['total_nts']
report['AtoC_freq'] = report['AtoC_nts'] / report['total_nts']
report['CtoA_freq'] = report['CtoA_nts'] / report['total_nts']
report['AtoT_freq'] = report['AtoT_nts'] / report['total_nts']
report['TtoA_freq'] = report['TtoA_nts'] / report['total_nts']
report['CtoG_freq'] = report['CtoG_nts'] / report['total_nts']
report['GtoC_freq'] = report['GtoC_nts'] / report['total_nts']
report['CtoT_freq'] = report['CtoT_nts'] / report['total_nts']
report['TtoC_freq'] = report['TtoC_nts'] / report['total_nts']
report['GtoT_freq'] = report['GtoT_nts'] / report['total_nts']
report['TtoG_freq'] = report['TtoG_nts'] / report['total_nts']
report.to_csv(od + exp + "_variant_summary.txt", sep="\t", index=False)
