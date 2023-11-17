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

from matplotlib.cm import get_cmap
from sklearn import preprocessing

import matplotlib.ticker as ticker


def get_information_model(model_result, target=None, analysis=None, gender=None, formula=None, alpha=0.01):
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
    result = result[['effect', 'odds_ratio', 'odds_ratio_lower', 'odds_ratio_upper', 'p_value']]

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

def run_regressions(dataset, targets, families, analysis_tag, genders, s_het_list=['s_het']):
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

        print (f"\tProcessing {gender} samples", flush=True)

        # select samples
        dataset_subset = select_samples_gender(dataset=dataset, gender=gender)

        for target, family in zip(targets, families):

            print (f"\t\tProcessing {target}", flush=True)
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
                                           formula=formula)


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
        ['effect', 'formula', 'odds_ratio_lower', 'odds_ratio_upper'], axis=1)

    # special case for edu-all
    df = df[~df['analysis'].str.contains("any education = 'all'")]

    # fix output format
    df["odds_ratio"] = df["odds_ratio"].apply(lambda x: round(x, 2))
    df["p_value_raw"] = df["p_value"].copy()
    df["p_value"] = df["p_value"].apply(sci_notation)
    df['n_observations'] = df['n_observations'].apply(lambda x: "{0:,}".format(x))

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

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

# Plots for paper
def get_plot_data(df, target, tag):
    """
        Selects regression results for a `target`, 
        converts odds ratio confidence intervals into distances to odds ratio (needed for plotting)
    """
    # select necessary target
    data = df[df['target'] == target].copy()

    # convert confidence intervals into distances
    data['odds_ratio_lower'] = data['odds_ratio'] - data['odds_ratio_lower']
    data['odds_ratio_upper'] = data['odds_ratio_upper'] - data['odds_ratio']

    # ad tag
    data['tag'] = tag

    return data[['target', 'odds_ratio', 'odds_ratio_lower', 'odds_ratio_upper', 'p_value_pretty', 'tag', 'gender', 'feature', 'analysis']]

def order_df(df, column, order):
    """
    Creates dataframe with values in `column` organized in `order`. 
    """
    df_new = []

    for value in order:
         df_new.append(df[df[column] == value])

    df_new = pd.concat(df_new)

    return df_new


def plot_errorbar_grouped(df, axis, y_column, group_column, title, 
                          xlim=None, ymargin=0.4, legend_loc='lower right', legend_kwargs={}, group_scale=0.1, y_scale=1,
                          y_order=None, group_order=None, vertical_loc=1, colors=get_cmap("Accent").colors):
    """
        Makes a horizontal errorbar plot for odds ratio, grouped by `group_column`.  
    """

    # we will use label encoder to calculate the positions
    # of the data on the plot
    le = preprocessing.LabelEncoder()

    # calculate positions of the "groups"
    if y_order is None:
        df['y_group_loc'] = le.fit_transform(df[y_column])*y_scale
    else:
        df['y_group_loc'] = df[y_column].apply(lambda x: y_order[::-1].index(x))*y_scale

    # calculate positions of each plot within "groups"
    if group_order is None:
        df['y'] = df['y_group_loc'] + le.fit_transform(df[group_column])*group_scale
    else:
        df['y'] = df['y_group_loc'] + df[group_column].apply(lambda x: group_order[::-1].index(x))*group_scale

    # plot vertical line, denoting absence of effect   
    axis.axvline(x=vertical_loc,  color='gray', linewidth=1, alpha=0.5, linestyle='dotted',)

    # calculate the positions of labels on y-axis
    y_labels_dict = df[[y_column, 'y']].groupby(y_column).mean().apply(lambda x: round(x, 2)).to_dict()['y']

    # sort values in our df by position, so we control
    # for the order in the legend
    df = df.sort_values(by=['y'], ascending=False)

    # plot each group errorbar separatly with its own color
    for idx, group in enumerate(df[group_column].drop_duplicates().values):

        df_group = df[df[group_column] == group]

        axis.errorbar(x=df_group['odds_ratio'], 
                      y=df_group['y'], 
                      xerr=df_group[['odds_ratio_lower', 'odds_ratio_upper']].values.T, 
                      label=group,
                      capsize=0, marker='s', markersize=3, markeredgecolor=colors[idx], 
                      markerfacecolor=colors[idx], linestyle='', elinewidth=1, color=colors[idx])


    # add labels on y-axis
    axis.set_yticks(list(y_labels_dict.values()))
    axis.set_yticklabels(list(y_labels_dict.keys()))

    # add margins on top and bottom
    axis.margins(y=ymargin)

    # Customize spines
    axis.spines['left'].set_color('black')
    axis.spines['bottom'].set_color('black')
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    # Add ticks
    axis.yaxis.set_ticks_position('left')
    axis.xaxis.set_ticks_position('bottom')
    axis.tick_params(which='major', width=1.00, length=2.5)
    axis.tick_params(which='minor', width=0.75, length=1.25)

    # set limit on x-axis
    if xlim is not None:
        axis.set_xlim(xlim)

    # set title and grid style
    axis.set_title(title)
    axis.grid(linestyle='dotted')   

    # set ticks and labels on x-axis
    axis.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    axis.set_xlabel('OR (99% CI)')

    # add legend
    axis.legend(loc=legend_loc, numpoints=1, **legend_kwargs)

    return df


def plot_errorbar_grouped_transposed(df, axis, y_column, group_column, title, 
                          xlim=None, ymargin=0.4, legend_loc='lower right', group_scale=0.1, y_scale=1,
                          y_order=None, group_order=None, vertical_loc=1, colors=get_cmap("Accent").colors):
    
    """
        Makes a vertical errorbar plot for odds ratio, grouped by `group_column`.  
    """
    # we will use label encoder to calculate the positions
    # of the data on the plot
    le = preprocessing.LabelEncoder()

    # calculate positions of the "groups"
    if y_order is None:
        df['y_group_loc'] = le.fit_transform(df[y_column])
    else:
        df['y_group_loc'] = df[y_column].apply(lambda x: y_order[::-1].index(x))

    # calculate positions of each plot within "groups"
    if group_order is None:
        df['y'] = df['y_group_loc']*y_scale + le.fit_transform(df[group_column])*group_scale
    else:
        df['y'] = df['y_group_loc']*y_scale + df[group_column].apply(lambda x: group_order[::-1].index(x))*group_scale

    # calculate the positions of labels on x-axis
    y_labels_dict = df[[y_column, 'y']].groupby(y_column).mean().apply(lambda x: round(x, 2)).to_dict()['y']

    # sort values in our df by position, so we control
    # for the order in the legend
    df = df.sort_values(by=['y'], ascending=True)

    # plot each group errorbar separatly with its own color
    for idx, group in enumerate(df[group_column].drop_duplicates().values):

        df_group = df[df[group_column] == group]

        axis.errorbar(y=df_group['odds_ratio'], 
                      x=df_group['y'], 
                      yerr=df_group[['odds_ratio_lower', 'odds_ratio_upper']].values.T, 
                      label=group,
                      capsize=0, marker='s', markersize=10, markeredgecolor=colors[idx], 
                      markerfacecolor=colors[idx], linestyle='', elinewidth=3, color=colors[idx])

    # plot vertical line, denoting absence of effect 
    axis.axhline(y=vertical_loc,  linestyle='--', color='salmon', linewidth=3)

    # add labels on x-axis
    axis.set_xticks(list(y_labels_dict.values()))
    axis.set_xticklabels(list(y_labels_dict.keys()))

    # add margins in the left and right 
    axis.margins(x=ymargin)

    # set limits on x-axis
    if xlim is not None:
        axis.set_ylim(xlim)

    # set title and grid style
    axis.set_title(title)
    axis.grid(linestyle='dotted')   

    # set ticks and labels on y-axis
    axis.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    axis.set_ylabel('OR ()', fontsize=20)

    # add legend
    axis.legend(loc=legend_loc, numpoints=1)
