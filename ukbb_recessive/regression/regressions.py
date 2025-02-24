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
import numpy as np

import statsmodels.api as sm
import statsmodels.formula.api as smf


def get_information_model(model_result, target=None, analysis=None, gender=None, formula=None, alpha=0.01, family=None, n_tests_correction=1):
    """
        Extracts effect estimation and p-values from fitted regression model, 
        calculates odds ratios and confidence intervals, annotates result with tags 
        (target, gender, analysis, formula)
    """
    # extract effects with confidence intervals and p_values
    # index conains information about the regressors
    result = pd.concat([model_result.params, model_result.conf_int(alpha=alpha), model_result.pvalues], axis=1)
    result.columns = ['effect', 'lower', 'upper', 'p_value']

    # remove birth cohort since np.exp goes to infinity
    result = result.iloc[~result.index.str.contains('birth_cohort'), :]
    
    # convert effects into odds ratios
    result['odds_ratio'] = np.exp(result['effect'])
    result['odds_ratio_lower'] = np.exp(result['lower'])
    result['odds_ratio_upper'] = np.exp(result['upper'])

    if (family == 'binomial'):
        result['odds_ratio_pretty'] = 'OR = ' + result['odds_ratio'].apply(lambda x: f"{x:.3f}")
        result[f'{100-alpha*100}% CI'] = result.apply(lambda x: f"[{x['odds_ratio_lower']:.3f}, {x['odds_ratio_upper']:.3f}]", axis=1)
    else: 
        result['odds_ratio_pretty'] = 'ES = ' + result['effect'].apply(lambda x: f"{x:.3f}") 
        result[f'{100-alpha*100}% CI'] = result.apply(lambda x: f"[{x['lower']:.3f}, {x['upper']:.3f}]", axis=1)

    result.loc[:, "p_value_corrected"] = result["p_value"]*n_tests_correction
    result.loc[:, 'bonferroni_correction_coef'] = n_tests_correction


    result = result[['effect', 'odds_ratio', 'odds_ratio_lower', 'odds_ratio_upper', 'odds_ratio_pretty', 
                     f'{100-alpha*100}% CI', 'p_value', 'p_value_corrected', 'bonferroni_correction_coef']]

    # add information about observations used for regression model
    result['n_observations'] = model_result.nobs
    
    # tag results
    if target is not None:
        result['target'] = target

    if gender is not None:
        result['gender'] = gender

    if analysis is not None:
        result['analysis'] = analysis
    
    if formula is not None:
        result['formula'] = formula

    if family is not None:
        result['family'] = family
    
    # add regressors as a column and return
        
    return result.reset_index().rename(columns={'index': 'feature'})


def get_formula(target, s_het_list):
    """
        Returns regression formula, can take several s_hets as an input
    """
    pca_columns = [f"PCA_{i}" for i in range (1, 41)]

    formula = f"{target} ~ {' + '.join(s_het_list)} + age_at_recruitment + I(age_at_recruitment**2)"
    formula = formula + ' + ' + ' + '.join(pca_columns) 

    return formula


def select_samples_gender(dataset, gender):
    """
        Filters dataset based on gender if gender in {'males', 'females'};
        doesn't apply any filter if gender = 'all',
    """
    if gender == 'all':
        return dataset
    elif gender == 'males':
        return dataset[dataset['gender'] == 1]
    elif gender == 'females':
        return dataset[dataset['gender'] == 0]
    else:
        raise Exception("sample_tag should be one of {'all', 'males', 'females'}")
    

def get_target_family(family):
    """
        Return corresponfing GLM family class based on string input. 
    """
    if family == 'binomial':
        return sm.families.Binomial()
    elif family == 'poisson':
        return sm.families.Poisson()
    elif family == 'gaussian':
        return sm.families.Gaussian()
    else:
        raise Exception("family should be one of {'binomial', 'poisson', 'gaussian'}")

def run_regressions(dataset, targets, families, analysis_tag, genders, s_het_list=['s_het'], tab_offset='', n_tests_correction=1):
    """
        Main function that runs regressions on `dataset`. 
        The function runs regression for all `targets` , `genders` one-by-one. 

        `targets` is a list of all phenotypes of interest. that would be considered as a response variable. 
        `families` is a list of all corresponding families used for a target in GLM. 
        `s_het_list` is a list of all regressors, that should be included into the model (except PCA and age that are inculded by default)
        `genders` is a list of all gender-specific filters; ['all'] will correspond to using unfiltedred dataset. 
        `analysis_tag` is a information, that would be added to the output of the model to identify the analysis.

        Returns pandas DataFrame with all regression results. 

    """

    assert len(targets) == len(families), 'Please check, that len(targets) == len(families)'

    # initialize results  
    results = []

    for gender in genders:

        print (f"{tab_offset}\tProcessing {gender} samples", flush=True)

        # select samples
        dataset_subset = select_samples_gender(dataset=dataset, gender=gender)

        for target, family in zip(targets, families):

            print (f"{tab_offset}\t\tProcessing {target}", flush=True)
            # generate formula
            formula = get_formula(target=target, s_het_list=s_het_list)

            # run regressions
            model = smf.glm(formula = formula, data=dataset_subset, family=get_target_family(family))
            fitted_model = model.fit()

            # getting effetcs + p_values
            result = get_information_model(fitted_model,
                                           target = target, 
                                           gender = gender, 
                                           analysis=analysis_tag, 
                                           formula=formula, 
                                           family=family, 
                                           n_tests_correction=n_tests_correction)


            # select interesting features
            result = result [result['feature'].apply(lambda x: ('s_het' in x) or (x == 'age_at_recruitment') or (x in s_het_list))]

            # save results
            results.append(result)

    results = pd.concat(results)

    return results

# pretty output formats

def SuperScriptinate(number):
  """
    Changes all digits to its superscript.
  """
  return number.replace('0','⁰').replace('1','¹').replace('2','²').replace('3','³').replace('4','⁴').replace('5','⁵').replace('6','⁶').replace('7','⁷').replace('8','⁸').replace('9','⁹').replace('-','⁻')

def NormalScriptinate(number):
  """
    Changes all superscript digits to its normal formart.
  """
  return number.replace('⁰', '0').replace('¹', '1').replace('²', '2').replace('³', '3').replace('⁴', '4').replace('⁵', '5').replace('⁶', '6').replace('⁷', '7').replace('⁸', '8').replace('⁹', '9').replace('⁻', '-')


def sci_notation(number, sig_fig=2):
  """
    Converts the number into the pretty text in scientific notation:
    For example: 1e-17 -> 1×10¹⁷
  """

  # convert numbet into string in scientific notation
  ret_string = "{0:.{1:d}e}".format(number, sig_fig)
  # divide into main and power: in 1e-17 1 is main, -17 is power
  main, power = ret_string.split("e")
  # removes leading "+" and strips leading zeros in power
  power = int(power)         
  # concatenate main and power into pretty string
  return main + "×10" + SuperScriptinate(str(power))

def p_value_sci_notation(number, sig_fig=2, treshold=0.01):
    """
        Converts the number into the pretty text in scientific notation:
        For example: 1e-17 -> 1×10¹⁷
    """
    if number > 1:
        return 1
    elif number < treshold:
        return sci_notation(number, sig_fig)
    else:
        return f"{round(number, 3)}"

def revert_sci_notation(number):
    """
        Changes all superscript digits to its normal formart.
    """
    return NormalScriptinate(number)


def prettify_table_for_paper(df, keep_effects=[]):
    """
        Creates table with relevant columns, rounded numbers and pretty scientific format for p-values
    """
    # drop unused columns
    df = df[df['feature'].str.contains('s_het') | df['feature'].isin(keep_effects)].drop(
        ['effect', 'formula', 'odds_ratio', 'odds_ratio_lower', 'odds_ratio_upper', 'family'], axis=1)

    # special case for edu-all
    df = df[~df['analysis'].str.contains("any education = 'all'")].copy()

    # fix output format
    # df.loc[:, "odds_ratio"] = df["odds_ratio"].apply(lambda x: round(x, 2))
    df.loc[:, "p_value_pretty"] = df["p_value"].apply(p_value_sci_notation).astype(str)
    df.loc[:, "p_value_corrected_pretty"] = df["p_value_corrected"].apply(p_value_sci_notation).astype(str)
    df.loc[:, 'n_observations_pretty'] = df['n_observations'].apply(lambda x: "{0:,}".format(x)).astype(str)
    

    df.drop(['n_observations', 'p_value', 'p_value_corrected', 'bonferroni_correction_coef'], axis=1, inplace=True)

    return df



def save_table_for_paper(all_results, path, keep_effects=[]):
    """
        Saves dictionary of DataFrames with results as Excel table. 
        Every key, value pair in the dictionary of results will be added as columns. 

        Prettifies the output result and saves as Excel with two sheets: raw and pretty. 

    """
    all_results_df = pd.concat(all_results.values(), keys=all_results.keys(), axis=1)

    # prettify the results
    all_results_pretty = {k: prettify_table_for_paper(v, keep_effects=keep_effects) for k,v in all_results.items()}
    all_results_pretty_df = pd.concat(all_results_pretty.values(), keys=all_results.keys(), axis=1)

    # create one column with this info
    some_key = list(all_results.keys())[0]
    all_results_pretty_df[('Info', 'feature')] = all_results_pretty_df[(some_key, 'feature')].copy()
    all_results_pretty_df[('Info', 'target')] = all_results_pretty_df[(some_key, 'target')].copy()
    all_results_pretty_df[('Info', 'gender')] = all_results_pretty_df[(some_key, 'gender')].copy()
    all_results_pretty_df[('Info', 'analysis')] = all_results_pretty_df[(some_key, 'analysis')].copy()

    # delete the rest of the columns
    drop_columns = []
    for key in all_results:
        drop_columns += [(key, 'target'), (key, 'gender'), (key, 'analysis')]
    all_results_pretty_df = all_results_pretty_df.drop(drop_columns, axis=1)

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(path, engine="xlsxwriter")

    # Write each dataframe to a different worksheet.
    all_results_df.to_excel(writer, sheet_name="Raw data")
    all_results_pretty_df.to_excel(writer, sheet_name="Paper table")

    # Bug/weird behaviour of MultiIndex in ExcelWriter
    writer.sheets["Raw data"].set_row(2, None, None, {'hidden': True})
    writer.sheets["Paper table"].set_row(2, None, None, {'hidden': True})

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()
    

def read_results_excel(path, flatten_multiindex=True):
    # read raw table
    reader = pd.ExcelFile(path)

    all_results_df = pd.read_excel(reader, sheet_name="Raw data", header=[0, 1], skiprows=[2])
    all_results_df = all_results_df.drop(all_results_df.columns[0], axis=1)

    # prettify p-values
    new_columns = [(level0, 'p_value_pretty') for level0 in all_results_df.columns.get_level_values(level=0).unique()]
    all_results_df[new_columns] = all_results_df.loc[:, (slice(None), 'p_value')].map(p_value_sci_notation)
    
    if 'p_value_corrected' in all_results_df.columns.get_level_values(level=1):
        new_columns = [(level0, 'p_value_corrected_pretty') for level0 in all_results_df.columns.get_level_values(level=0).unique()]
        all_results_df[new_columns] = all_results_df.loc[:, (slice(None), 'p_value_corrected')].map(p_value_sci_notation)

    if flatten_multiindex:
        all_results_df = all_results_df.stack(level=0, future_stack=True).reset_index().drop('level_0', axis=1).rename({'level_1': 'dataset'}, axis=1)

    return all_results_df
