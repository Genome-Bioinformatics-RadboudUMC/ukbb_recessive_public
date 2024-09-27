'''
This file is part of the "ukbb_recessive" project, that explores the 
negative selection on heterozygotes and the landscape of recessive disease
using UKBB data. 

Copyright (C) 2023  Gelana Khazeeva
Copyright (C) 2023  Hila Fridman

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

from csv import DictReader, DictWriter
import tqdm
import pandas as pd
import numpy as np

class VariantFeatures:
    def __init__(self):
        pass

    def filter_plps_in_samples(self, rap_files, output_folder, all_plps_file):
        """
        This function leaves only PLPs in samples full variants list.
        :param rap_files: list of files to filter
        :param output_folder: output folder
        :param all_plps_file: path to all PLPs file
        :return: nothing
        """
        print()
        print ("Entering `filter_plps_in_samples` function...")

        # output fieldnames
        fieldnames = ['chrom', 'pos', 'ref', 'alt', 'GT', 's',
                      'raw pos', 'raw ref', 'raw alt', 'alleles', 'GT.alleles']

        # read all PLPs file
        plps_set = create_plps_set(all_plps_file)
        print (f"Number of total PLPs: {len(plps_set)}")

        # filter PLPs in samples files
        for rap_file in tqdm.tqdm(rap_files):

            output_file = f"{output_folder}/{get_file_name(rap_file)}"

            with open(rap_file, 'r') as f_in, open(output_file, 'w') as f_out:
                reader = DictReader(f_in)
                writer = DictWriter(f_out, fieldnames=fieldnames)
                writer.writeheader()

                for line in reader:
                    if (line['chrom'], line['pos'], line['ref'], line['alt']) in plps_set:
                        writer.writerow(line)

        print(f"Function `filter_plps_in_samples` finished, result written in `{output_folder}`.")
        print()

    def collect_rare_plps(self, het_occurence_threshold, all_plps_file, s_het_file, genes_list=None, hom_occurence_threshold=0.):
        """
        Creates DataFrame with rare PLPs according to the `het_occurrence_threshold`, adds information about s_het
        :param het_occurence_threshold: threshold for the heterozygous variant occurrence in the cohort
        :param all_plps_file: path to all PLPs file
        :param s_het_file: path to file with s_het scores for genes
        :param genes_list: list of genes, that need to be used
        :return: returns a DataFrame with rare PLPs annotated with s_het scores
        """
        print()
        print ("Entering `collect_rare_plps` function...")

        # read s_het data
        s_het_cassa = pd.read_csv(s_het_file,sep='\t').rename(columns={'gene_symbol': 'gene'}).drop_duplicates()
        assert s_het_cassa.shape[0] == s_het_cassa['gene'].drop_duplicates().shape[0], "S-het file contains duplicates!"

        # read all plps
        all_plps = pd.read_csv(all_plps_file, sep='\t')
        print(f"Initial total numbers of PLPs: {all_plps.shape[0]}")

        # select rare PLPS
        rare_plps = all_plps[(all_plps['hets'] <= het_occurence_threshold) & (all_plps['homs'] <= hom_occurence_threshold)]
        rare_plps = rare_plps[['chr', 'position', 'ref',
                               'alt', 'gene', 'hets', 'homs']].rename(columns={'chr': 'chrom', 'position': 'pos'})
        print(f"Total numbers rare PLPs using <treshold={het_occurence_threshold}>: {rare_plps.shape[0]}")

        # add s_het informations
        rare_plps = rare_plps.merge(s_het_cassa, how='left', on='gene')

        # filter_out gene panel if necessary
        if genes_list is not None:
            genes_list = set(genes_list)
            rare_plps = rare_plps[rare_plps['gene'].astype(str).isin(genes_list)]
            print(f"Total numbers rare PLPs in specified gene list: {rare_plps.shape[0]}")

        print(f"Function `collect_rare_plps` finished.")
        print()

        return rare_plps

    def read_sample_plps(self, cohort_plps_files, filter_homozygous=True):
        """
        reads samples PLPs and merges them together
        :param cohort_plps_files: list of all plp files for cohort
        :return: DataFrame with all data merged together
        """
        plps = []
        for plp_file in cohort_plps_files:
            try:
                plps.append(pd.read_csv(plp_file))
            except:
                continue

        plps = pd.concat(plps).drop_duplicates()
        print ("All PLPs in the cohort:", plps.shape[0])
        
        if filter_homozygous:
            # remove homozygous carriers
            plps = plps[plps['GT'] != '1/1']
            print ("Heterozygous PLPs in the cohort:", plps.shape[0])

        return plps

    def calculate_s_het_per_sample(self, het_occurence_threshold, all_plps_file, s_het_file,
                                   cohort_plps_files, genes_list=None, high_shet_threshold=0.15):
        """
        Calculates s_het and number of mutations per sample
        :param het_occurence_threshold: threshold for the heterozygous variant occurrence in the cohort
        :param all_plps_file: path to all PLPs file
        :param s_het_file: path to file with s_het scores for genes
        :param cohort_plps_files: path to PLPs for the cohort
        :return: DataFrame with s_het score and mutation counts per sample
        """
        print()
        print ("Entering `calculate_s_het_per_sample` function...")

        # select rare PLPs
        rare_plps = self.collect_rare_plps(het_occurence_threshold=het_occurence_threshold,
                                           all_plps_file=all_plps_file,
                                           s_het_file=s_het_file,
                                           genes_list=genes_list)

        # read cohort PLPs
        plps = self.read_sample_plps(cohort_plps_files)
        print(f"Total numbers of PLP variants in cohort: {plps.shape[0]}")

        # filter rare PLPs in cohort
        plps = plps.merge(rare_plps)
        print(f"Total numbers of rare PLP variants in cohort: {plps.shape[0]}")
        
        #add a column for high shet gene yes/no
        plps['high_shet'] = (plps['s_het'] >= high_shet_threshold).astype(int)
        plps['high_shet'] = plps['high_shet'].fillna(0)

        # add a column for a non-additive s_het
        plps['max_s_het'] = plps['s_het'].copy()

        # several PLPs in the same gene are treated as one
        plps = plps.groupby(['s', 'gene'])[['chrom', 's_het', 'max_s_het', 'high_shet']].max().reset_index()
        print(f"Total numbers of rare PLP  variants in cohort, one per gene: {plps.shape[0]}")

        # calculate s_het and mutation counts per sample
        plps_shet = plps.groupby('s').agg({'s_het': get_cumm_s_het,
                                           'chrom': 'count',
                                           'max_s_het': 'max',
                                           'high_shet': 'max'
                                           }).reset_index().rename(columns={'s': 'eid', 
                                                                            'chrom':'mutations_cnt', 
                                                                            'high_shet':'has_high_shet_var'})

        plps_shet['has_mutation'] = (plps_shet['mutations_cnt'] > 0).astype(int)

        print(f"Function `calculate_s_het_per_sample` finished.")
        print()

        return plps_shet
    


class VariantFeaturesLoF:
    def __init__(self):
        pass

    def filter_lofs_in_samples(self, rap_files, output_folder, all_lof_file):
        """
        This function leaves only LoF in samples full variants list.
        :param rap_files: list of files to filter
        :param output_folder: output folder
        :param all_lof_file: path to all singletons file
        :return: nothing
        """
        print()
        print ("Entering `filter_lofs_in_samples` function...")

        # output fieldnames
        fieldnames = ['s', 'chrom', 'pos', 'ref', 'alt', 'updated GT','locus','alleles','GT']

        # read all PLPs file
        lofs_set = create_lofs_set(all_lof_file)
        print (f"Number of total LoFs: {len(lofs_set)}")

        # filter LoFs in samples files
        for rap_file in tqdm.tqdm(rap_files):

            output_file = f"{output_folder}/{get_file_name(rap_file)}"

            with open(rap_file, 'r') as f_in, open(output_file, 'w') as f_out:
                reader = DictReader(f_in)
                writer = DictWriter(f_out, delimiter='\t', fieldnames=fieldnames)
                writer.writeheader()

                for line in reader:
                    if (line['chrom'], line['pos'], line['ref'], line['alt']) in lofs_set:
                        writer.writerow(line)

        print(f"Function `filter_lofs_in_samples` finished, result written in `{output_folder}`.")
        print()

    def collect_lofs(self, all_lofs_file, s_het_file, genes_list=None):
        """
        Creates DataFrame with LoFs from all or necessary genes if 'gene_list' provided, adds information about s_het
        :param all_lofs_file: path to all LoFs file
        :param s_het_file: path to file with s_het scores for genes
        :param genes_list: list of genes, that need to be used
        :return: returns a DataFrame with LoFs annotated with s_het scores
        """
        print()
        print ("Entering `collect_lofs` function...")

        # read s_het data
        s_het = pd.read_csv(s_het_file, sep='\t').rename(columns={'gene_symbol': 'gene'}).drop_duplicates()

        # read all plps
        all_lofs = pd.read_csv(all_lofs_file, sep='\t').rename(columns={'SYMBOL': 'gene'})
        all_lofs[['ref', 'alt']] =  all_lofs[['ref', 'alt']].fillna('')
        print(f"Initial total numbers of LoFs: {all_lofs.shape[0]}")

        # add s_het information
        all_lofs = all_lofs.merge(s_het, how='left', on='gene')

        # filter_out gene panel if necessary
        if genes_list is not None:
            genes_list = set(genes_list)
            all_lofs = all_lofs[all_lofs['gene'].astype(str).isin(genes_list)]
            print(f"Total numbers rare LoF in specified gene list: {all_lofs.shape[0]}")


        # handling gene overlap (same variant could be annotated with more than 1 gene)
        all_lofs = (
            all_lofs
            .sort_values(by=['chrom', 'pos', 'ref', 'alt', 's_het', 'gene'], ascending=False)
            .groupby(['chrom', 'pos', 'ref', 'alt'])
            .agg({'s_het': 'max', 'gene': 'first'})
            .reset_index()
        )

        print(f"Total numbers of LoFs after handling gene overlap: {all_lofs.shape[0]}")
        print(f"Total numbers of LoFs with null s_het after handling gene overlap: {all_lofs['s_het'].isnull().sum()}")
        print(f"Function `collect_lofs` finished.")
        print()

        return all_lofs

    def read_sample_lofs(self, cohort_lofs_files):
        """
        reads samples LoFs and merges them together
        :param cohort_lofs_files: list of all lofs files for cohort
        :return: DataFrame with all data merged together
        """
        lofs = []
        for lof_file in cohort_lofs_files:
            try:
                lofs.append(pd.read_csv(lof_file, sep='\t'))
            except:
                continue

        lofs = pd.concat(lofs).drop_duplicates()
        lofs[['ref', 'alt']] =  lofs[['ref', 'alt']].fillna('')
        print ("All LOFs in the cohort:", lofs.shape[0])

        # remove homozygous carriers
        lofs = lofs[lofs['updated GT'] != '1/1']
        print ("Heterozygous-only LOFs in the cohort:", lofs.shape[0])

        return lofs

    def calculate_s_het_per_sample(self, all_lofs_file, s_het_file,
                                   cohort_lofs_files, genes_list=None, high_shet_threshold=0.15):
        """
        Calculates s_het and number of mutations per sample
        :param all_lofs_file: path to all LoF file
        :param s_het_file: path to file with s_het scores for genes
        :param cohort_lofs_files: path to LoFs for the cohort
        :return: DataFrame with s_het score and mutation counts per sample
        """
        print()
        print ("Entering `calculate_s_het_per_sample` function...")

        # select LoFs
        lofs = self.collect_lofs(all_lofs_file=all_lofs_file, 
                                 s_het_file=s_het_file, 
                                 genes_list=genes_list)

        # read cohort LoFs
        cohort_lofs = self.read_sample_lofs(cohort_lofs_files)
        print(f"Total numbers of LoF variants in cohort: {cohort_lofs.shape[0]}")

        # filter LoFs in cohort
        cohort_lofs = cohort_lofs.merge(lofs,  on=['chrom', 'pos', 'ref', 'alt']).drop_duplicates()
        print(f"Total numbers of selected LoF variants in cohort: {cohort_lofs.shape[0]}")
        
        #add a column for high shet gene yes/no
        cohort_lofs['high_shet'] = (cohort_lofs['s_het'] >= high_shet_threshold).astype(int)
        cohort_lofs['high_shet'] = cohort_lofs['high_shet'].fillna(0)

        # add a column for a non-additive s_het
        cohort_lofs['max_s_het'] = cohort_lofs['s_het'].copy()

        # several LOFs in the same gene are treated as one
        cohort_lofs = cohort_lofs.groupby(['s', 'gene'])[['chrom', 's_het', 'max_s_het', 'high_shet']].max().reset_index()
        print(f"Total numbers of selected LoF variants in cohort, one per gene: {cohort_lofs.shape[0]}")

        # calculate s_het and mutation counts per sample
        lofs_shet = cohort_lofs.groupby('s').agg({'s_het': get_cumm_s_het,
                                                  'chrom': 'count',
                                                  'max_s_het': 'max',
                                                  'high_shet': 'max'
                                                  }).reset_index().rename(columns={'s': 'eid', 
                                                                                   'chrom':'mutations_cnt', 
                                                                                   'high_shet':'has_high_shet_var'})

        lofs_shet['has_mutation'] = (lofs_shet['mutations_cnt'] > 0).astype(int)

        print(f"Function `calculate_s_het_per_sample` finished.")
        print()

        return lofs_shet

def get_cumm_s_het(s_het_list):
    """
        calculate s_het for a sample based on PLPs
    """
    return 1 - np.prod([1-np.array(s_het_list)])

def get_file_name(path):
    return path.split('/')[-1]

def create_plps_set(plps_file):
    """
        read PLPs using DictReader and converts this list to set
    """
    plps = []

    with open(plps_file, 'r') as f_in:
        reader = DictReader(f_in, delimiter='\t')

        for line in reader:
            plps.append(tuple([line['chr'], line['position'], line['ref'], line['alt']]))

    plps = set(plps)

    return plps


def create_lofs_set(lofs_file):
    """
        read LoF variants using DictReader and converts this list to set
    """
    lofs = []

    with open(lofs_file, 'r') as f_in:
        reader = DictReader(f_in, delimiter='\t')

        for line in reader:
            lofs.append(tuple([line['chrom'], line['pos'], line['ref'], line['alt']]))

    lofs = set(lofs)

    return lofs
