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

from ukbb_recessive.data_collection.phenotypes import PhenotypeFeatures
from ukbb_recessive.data_collection.variants import VariantFeatures, VariantFeaturesLoF

class RegressionDataset:
    """
        This class contains functions for the dataset collection for 
        a regression analysis. 
    """
    def __init__(self, age_children_path=None, pca_path=None, other_features_path=None,
                 het_occurrence_threshold=None, all_plps_file=None, s_het_file=None,
                 cohort_plps_files=None, samples_list=None, genes_list=None, 
                 high_shet_threshold=0.15, dataset='recessive'):
        
        # phenotype features paths (from the UK Biobank)
        self.age_children_path = age_children_path
        self.pca_path = pca_path
        self.other_features_path = other_features_path

        # variant features params
        self.het_occurrence_threshold = het_occurrence_threshold # maximum heterozyguos occurrence in the cohort  
        self.all_plps_file = all_plps_file # path to a file with all variants of interest 
        self.s_het_file = s_het_file # path to a file with a selective constraint
        self.cohort_plps_files = cohort_plps_files # paths to files with cohort variants information
        self.high_shet_threshold = high_shet_threshold # redundant, but this is a treshold to consider a gene as highly constraint
        self.dataset = dataset # could be `synonymous`, `recessive` or `lof`, defines which variant features class to invoke (PLPs or LOFs)

        #filters
        self.samples_list = samples_list # list of samples of interest
        self.genes_list = genes_list #list of genes of interest (used to calculate s-het over them)

    def collect_phenotypic_features(self):
        """
            This function collects phenotypic features, based on UKBB data. 
            It will leave only samples that are defined in `self.samples_list`.
        """
        # collect phenotypic features
        features = PhenotypeFeatures().collect_phenotype_features(
            age_children_path = self.age_children_path,
            pca_path = self.pca_path,
            other_features_path = self.other_features_path
        )

        # filter-out samples
        if self.samples_list is not None:
            features = features[features['eid'].astype(str).isin(self.samples_list)]
            print(f"Phenotypic features after filtration: {features.shape}")

        # rename feature as glm doesnt allow "/" in names
        features = features.rename(columns={'uni_1/0_including_none':'uni_1_0_including_none'})

        return features
    
    def collect_variant_features(self):
        """
            This function collects variant-based features, which is:
            genetic (s-het) burden, number of affected genes, carriership status etc.
        """
        print ("Dataset type:", self.dataset)

        # this code is used for PLPs in recessive genes
        if (self.dataset == 'recessive') or (self.dataset == 'synonymous'):
            # collect variant features
            variant_features = VariantFeatures().calculate_s_het_per_sample(
                het_occurence_threshold = self.het_occurrence_threshold,
                all_plps_file = self.all_plps_file,
                s_het_file = self.s_het_file,
                cohort_plps_files = self.cohort_plps_files,
                genes_list = self.genes_list,
                high_shet_threshold = self.high_shet_threshold
                
            )
        # this code is used for singleton LoFs
        elif self.dataset == 'lof':
            print (f"Ignoring het_occurence_threshold = {self.het_occurrence_threshold}")
            
            variant_features = VariantFeaturesLoF().calculate_s_het_per_sample(
                all_lofs_file = self.all_plps_file, 
                s_het_file=self.s_het_file,
                cohort_lofs_files=self.cohort_plps_files, 
                genes_list=self.genes_list, 
                high_shet_threshold=self.high_shet_threshold
            )

        else:
            raise Exception("dataset should be one of ('synonymous', 'recessive', 'lof')")
        
        return variant_features


    def collect(self):
        """
            This function is redundant but it allows to create joint dataset
            of phenotypic and variant features together. 
        """
        print()
        print("Entering `RegressionDataset.collect` function...")

        features = self.collect_phenotypic_features()
        variant_features = self.collect_variant_features()

        # merge features together
        dataset = features.merge(variant_features, on='eid', how='left')
        print(f"Total dataset for regression shape: {dataset.shape}")

        print(f"Function `RegressionDataset.collect` finished.")
        print()

        return dataset
