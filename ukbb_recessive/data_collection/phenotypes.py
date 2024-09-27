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

import pandas as pd
from math import nan
import numpy as np

class PhenotypeFeatures:
    """
    Incapsulates the phenotypic feature collection process for UKBB
    """
    def __init__(self):
        self.features_meaning = {
            '34-0.0': 'birth_cohort',
            '31-0.0': 'gender',
            '2405-0.0': 'number_of_children_fathered',
            '2734-0.0': 'number_of_live_births',
            '21022-0.0': 'age_at_recruitment',
            '6138-0.0': 'qualifications',
            '709-0.0': 'people_number_in_household',
            '6141-0.0': 'people_related_in_household',
            '41202-0.0': 'diagnosis_main_ICD10',
            '41204-0.0': 'diagnosis_secondary_ICD10',
            '20544-0.0': 'mental_health_problems',
            '20005-0.0': 'email',
            '22006-0.0': 'genetic_ethnic_group',
            '20016-0.0': 'fluid_intelligence_score',
            '2139-0.0': 'age_first_sex',
            '42040-0.0': 'gp_record',
            '1707-0.0': 'handedness',
            '1747-0.0': 'hair_color',
            '50-0.0': 'height',
            'p26410' : 'multiple_deprivation_engand',
            'p26427' : 'multiple_deprivation_scotland',
            'p26426' : 'multiple_deprivation_wales',
            'p26414' : 'edu_deprivation_england',
            'p26431' : 'edu_deprivation_scotland',
            'p26421' : 'edu_deprivation_wales',
            'p26411' : 'income_deprivation_england',
            'p26428' : 'income_deprivation_scotland',
            'p26418' : 'income_deprivation_wales',
            'p26415' : 'housing_deprivation_england',
            'p26432' : 'housing_deprivation_scotland',
            'p26423' : 'housing_deprivation_wales',
            'p26413' : 'health_deprivation_england',
            'p26430' : 'health_deprivation_scotland',
            'p26420' : 'health_deprivation_wales',
            'p738_i0' : 'household_income',
            'p1309_i0' : 'fresh_fruit_intake'
        }

    def _process_pca_features(self, pca_features):
        """
        converts pca feature names into PCA_1, PCA_2, ....
        :param pca_features: DataFrame with raw features
        :return: DataFrame with new feature names
        """

        pca_features.columns = [(f"PCA_{x.split('.')[-1]}" if '22009' in x else x) for x in pca_features.columns]

        return pca_features

    def collect_phenotype_features(self, age_children_path, pca_path, other_features_path):
        """
        collects all phenotypic features for the regression
        :param age_children_path: path to features with age and children info
        :param pca_path: path to the features with PCA
        :param other_features_path: path to the rest of the features
        :return:  DataFrame with all the features for every sample
        """
        print()
        print("Entering `collect_phenotype_features` function...")

        # read data
        age_children = pd.read_csv(age_children_path, sep='\t')
        pca = pd.read_csv(pca_path, sep='\t')
        other = pd.read_csv(other_features_path, sep='\t')

        # process pca features
        pca = self._process_pca_features(pca_features=pca)

        # merge the data together
        features = age_children.merge(other, on='eid', how='outer').merge(pca, on='eid', how='outer')

        # rename columns
        features = features.rename(columns=self.features_meaning)

        # print infortmation
        print("Current columns names list:", features.columns.tolist())
        print()
        print("Number of samlples with features:", features.shape[0])

        # birth_cohort
        features['birth_cohort'] = features['birth_cohort'].fillna(-1).apply(calc_birth_cohort)

        # number of childrens
        features.loc[(features['number_of_children_fathered'] < 0) |
                     (features['number_of_children_fathered'] > 7), 'number_of_children_fathered'] = None

        features.loc[(features['number_of_live_births'] < 0) |
                     (features['number_of_live_births'] > 7), 'number_of_live_births'] = None

        features['number_of_children_MF'] = features[['number_of_children_fathered',
                                                      'number_of_live_births']].max(axis=1)

        # childlessness
        features['childlessness'] = 1 - (features['number_of_children_MF'] > 0).fillna(0).astype(int)
        features.loc[features['number_of_children_MF'].isnull(), 'childlessness'] = None

        # university degree
        features['qualifications'] = features['qualifications'].fillna('')
        features["years_of_edu"] = features['qualifications'].fillna('').apply(years_of_edu)
        features['uni_1/0_excluding_none'] = features['qualifications'].apply(has_university_degree)
        features['uni_1/0_including_none'] = features['qualifications'].apply(has_university_degree_include_none)
        features['higher_education_including_none'] = features['qualifications'].apply(has_higher_education_include_none)
        features['any_education_including_none'] = features['qualifications'].apply(has_any_education_include_none)

        # living with a partner
        features['living_with_a_partner'] = ((features['people_number_in_household'] > 1) &
                                             (features['people_related_in_household'].str.contains("1"))).astype(int)
        features.loc[(features['people_number_in_household'] < 0) | features['people_number_in_household'].isnull(),
                    'living_with_a_partner'] = None

        # ICD mental health
        features['ICD_mental_health_yes_no'] = (features['diagnosis_main_ICD10'].apply(has_ICD_mental_health) |
                                                features['diagnosis_secondary_ICD10'].apply(has_ICD_mental_health))

        # mental health questionnarie
        features['mental_health_Q'] = (
            features['mental_health_problems']
                .apply(lambda x: str(x).split('|') if x is not None else [])
                .apply(lambda x: len(set(x) - {"-818", "-819"}) > 0).astype(int)
        )
        features.loc[features['mental_health_problems'].isnull(), 'mental_health_Q'] = None
        features.loc[features['mental_health_problems'].astype(str) == '', 'mental_health_Q'] = None

        # email
        features['email'] = features["email"].apply(lambda x: 0 if x == 2 else x)

        # ever_had_sex
        features['ever_had_sex'] = None
        features.loc[features['age_first_sex'] > 0., 'ever_had_sex'] = 1
        features.loc[features['age_first_sex'] == -2., 'ever_had_sex'] = 0

        # gp records
        features['has_gp_record'] = (features['gp_record'] > 0).astype(int)
        # features.loc[features['gp_record'].isnull(), 'has_gp_record'] = None

        # handedness
        features['is_left_handed'] = (features['handedness'] == 2.).astype(int)
        features.loc[features['handedness'].isin([3, -3]), 'is_left_handed'] = None

        # hair color
        features['is_blond'] = (features['hair_color'] == 1.).astype(int)
        features.loc[features['hair_color'] < 0, 'is_blond'] = None      

        # ICD diagnoses
        features['diagnosis_main_ICD10_cnt'] = features['diagnosis_main_ICD10'].apply(number_ICD_diagnoses)
        features['diagnosis_secondary_ICD10_cnt'] = features['diagnosis_secondary_ICD10'].apply(number_ICD_diagnoses)
        features['diagnosis_total_ICD10_cnt'] = features['diagnosis_main_ICD10_cnt'] + features['diagnosis_secondary_ICD10_cnt']

        features['diagnosis_secondary_ICD10_cnt_log'] = np.log(features['diagnosis_secondary_ICD10_cnt'])
        features['diagnosis_main_ICD10_cnt_log'] = np.log(features['diagnosis_main_ICD10_cnt'])
        features['diagnosis_total_ICD10_cnt_log'] = np.log(features['diagnosis_total_ICD10_cnt'])

        # infertility
        features['ICD_infertility'] = (features['diagnosis_main_ICD10'].apply(has_infertility_ICD) |
                                       features['diagnosis_secondary_ICD10'].apply(has_infertility_ICD)).astype(int)


        print(f"Function `collect_phenotype_features` finished.")
        print()

        return features


def calc_birth_cohort(year):
    """
    This function defines the birth cohort.
    """
    if year == -1:
        return None

    middle = int(year // 10 * 10 + 4)

    if year <= middle:
        return f"{middle - 4}-{middle}"
    else:
        return f"{middle + 1}-{middle + 5}"

def years_of_edu(qualifications):
    # based on:
    # https://academic.oup.com/ije/article/51/3/885/6521336#369530775
    qualifications = qualifications.split('|')
    
    if '1' in qualifications:
        return 20
    elif '5' in qualifications:
        return 19
    elif '6' in qualifications:
        return 15
    elif '2' in qualifications:
        return 13
    elif ('3' in qualifications) or ('4' in qualifications):
        return 10
    elif '-7' in qualifications:
        return 7
    else:
        return None

def has_university_degree(qualification):
    qualification = str(qualification).strip()
    
    # prefer not to answer or none of the above or empty
    if (qualification == '-3') | (qualification == '-7') | (qualification == ''):
        return None

    return int("1" in qualification.split("|"))


def has_university_degree_include_none(qualification):
    qualification = str(qualification).strip()

    # prefer not to answer or empty
    if (qualification == '-3')  | (qualification == ''):
        return None

    return int("1" in qualification.split("|"))


def has_higher_education_include_none(qualification):
    qualification = str(qualification).strip()

    # prefer not to answer or empty
    if (qualification == '-3')  | (qualification == ''):
        return None
    
    education = ["1","2","4","5","6"]
    
    for i in qualification.split("|"):
        if i in education:
            return 1

    return 0
        
    
def has_any_education_include_none(qualification):
    qualification = str(qualification).strip()

    # prefer not to answer or empty
    if (qualification == '-3')  | (qualification == ''):
        return None
    
    education = ["1","2","3","4","5","6"]
    
    for i in qualification.split("|"):
        if i in education:
            return 1
      
    return 0


def has_ICD_mental_health(codes):
    if (codes is None) or (str(codes) == ""):
        return None

    ICD_mental_health_codes = ["F200", "F201", "F202", "F203", "F204", "F205",
                               "F206", "F208", "F209", "F231", "F232", "F250",
                               "F251", "F252", "F258", "F259", "F840", "F841",
                               "F845", "F849", "F300", "F301", "F302", "F308",
                               "F309", "F310", "F311", "F312", "F313", "F314",
                               "F315", "F316", "F317", "F318", "F319", "F900",
                               "F909", "F700", "F701", "F708", "F709", "F711",
                               "F719", "F729", "F780", "F789", "F790", "F798",
                               "F799", "F800", "F801", "F802", "F803", "F809",
                               "F810", "F812", "F819", "F82", "F83", "F89"]

    ICD_mental_health_codes = set(ICD_mental_health_codes)

    codes_list = set(str(codes).split("|"))

    intersection_mental = codes_list.intersection(ICD_mental_health_codes)

    return len(intersection_mental) > 0


def number_ICD_diagnoses(codes):
    if (codes is None) or (str(codes) == ""):
        return 0
    codes_list = set(str(codes).split("|"))

    return len(codes_list)

def has_infertility_ICD(codes):
    infertility = set(['N46'] + [f"N97{i}" for i in range(10)])

    codes_list = set([code.split('.')[0] for code in str(codes).split("|")])

    intersection_infertility = codes_list.intersection(infertility)

    return len(intersection_infertility) > 0
        