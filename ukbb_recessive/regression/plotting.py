import matplotlib.pyplot as plt
import matplotlib
from matplotlib import font_manager

import pandas as pd
import numpy as np

from matplotlib.cm import get_cmap
from sklearn import preprocessing

import matplotlib.ticker as ticker

SMALL_SIZE = 5
MEDIUM_SIZE = 6
BIGGER_SIZE = 7

def add_fonts(font_dirs):
    font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

    for font_file in font_files:
        font_manager.fontManager.addfont(font_file)
        print ("Added:", font_file)

def configure_matplotlib():
    plt.rc('font', size=BIGGER_SIZE, family='Arimo') # controls default text sizes

    plt.rcParams['text.usetex']= False

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

def configure_axis(ax, ytick_size=MEDIUM_SIZE, xtick_size=MEDIUM_SIZE, 
                   xlabel_size=MEDIUM_SIZE, ylabel_size=MEDIUM_SIZE, 
                   x_label=None, y_label=None, ymargin=0.4, xmargin=None, 
                   xlim=None, ylim=None, 
                   format_x=False, format_y=False):
    
    # add margins on top and bottom
    if ymargin:
        ax.margins(y=ymargin)

    # add margins on left and right
    if xmargin:
        ax.margins(x=xmargin)

    # Turn off grid
    ax.grid(False) 

    # Customize spines
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

        
    # Set tick labels size
    ax.tick_params(axis='y', labelsize=ytick_size) 
    ax.tick_params(axis='x', labelsize=xtick_size) 

    # Add tick marks
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    ax.tick_params(which='major', width=1.00, length=2.5)
    ax.tick_params(which='minor', width=0.75, length=1.25)

    # add format
    if format_x:
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

    if format_y:
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

    # set ticks and labels on x-axis
    ax.set_xlabel(x_label, size=xlabel_size)
    ax.set_ylabel(y_label, size=ylabel_size)

    # set limit on x-axis
    if xlim is not None:
        ax.set_xlim(xlim)

    # set limit on x-axis
    if ylim is not None:
        ax.set_ylim(ylim)


# # Plots for paper
# def get_plot_data(df, target, tag):
#     """
#         Selects regression results for a `target`, 
#         converts odds ratio confidence intervals into distances to odds ratio (needed for plotting)
#     """
#     # select necessary target
#     data = df[df['target'] == target].copy()

#     # convert confidence intervals into distances
#     data.loc[:, 'odds_ratio_lower'] = data['odds_ratio'] - data['odds_ratio_lower']
#     data.loc[:, 'odds_ratio_upper'] = data['odds_ratio_upper'] - data['odds_ratio']

#     # ad tag
#     data.loc[:, 'tag'] = tag

#     return data[['target', 'odds_ratio', 'odds_ratio_lower', 'odds_ratio_upper', 'p_value_pretty', 'tag', 'gender', 'feature', 'analysis']]

def add_odds_ratio_intervals(data):
    """
        Converts odds ratio confidence intervals into distances to odds ratio (needed for plotting)
    """
    data.loc[:, 'odds_ratio_lower_distance'] = data['odds_ratio'] - data['odds_ratio_lower']
    data.loc[:, 'odds_ratio_upper_distance'] = data['odds_ratio_upper'] - data['odds_ratio']

    return data

def add_effect_size_intervals(data):
    """
        Converts odds ratio confidence intervals into distances to odds ratio (needed for plotting)
    """
    data.loc[:, 'effect_lower_distance'] = data['effect'] - np.log(data['odds_ratio_lower'])
    data.loc[:, 'effect_upper_distance'] = np.log(data['odds_ratio_upper']) - data['effect']

    return data

def order_df(df, column, order):
    """
    Creates dataframe with values in `column` organized in `order`. 
    """
    df_new = []

    for value in order:
         df_new.append(df[df[column] == value])

    df_new = pd.concat(df_new)

    return df_new

def encode_labels(column, order=None, scale=1):
    """
    Encodes labels in `column` to be used in plotting. 
    """
    le = preprocessing.LabelEncoder()

    # random order
    if order is None:
        return le.fit_transform(column)*scale
    
    # we want first value to be at the top
    return column.apply(lambda x: order[::-1].index(x))*scale


def plot_errorbar_grouped(df, axis, y_column, group_column,
                          legend_loc='lower right', legend_kwargs={}, group_scale=0.1, y_scale=1,
                          y_order=None, group_order=None, vertical_loc=1, colors=get_cmap("Accent").colors, 
                          plot_entity='odds_ratio', horizontal=False, 
                          markersize=3):
    """
        Makes a horizontal errorbar plot for odds ratio, grouped by `group_column`.  
    """

    # calculate positions of data points in the plot
    df.loc[:, 'y_group_loc'] = encode_labels(df[y_column], order=y_order, scale=y_scale)
    df.loc[:, 'y'] = df['y_group_loc'] + encode_labels(df[group_column], order=group_order, scale=group_scale)


    # plot vertical line, denoting absence of effect   
    if not horizontal:
        axis.axvline(x=vertical_loc, color='gray', linewidth=1, alpha=0.5, linestyle='dotted')
    else:
        axis.axhline(y=vertical_loc, color='gray', linewidth=1, alpha=0.5, linestyle='dotted')

    # calculate the positions of labels on y-axis
    y_labels_dict = df[[y_column, 'y']].groupby(y_column).mean().apply(lambda x: round(x, 2)).to_dict()['y']

    # sort values in our df by position, so we control
    # for the order in the legend
    df = df.sort_values(by=['y'], ascending=False)

    # plot each group errorbar separatly with its own color
    for idx, group in enumerate(df[group_column].drop_duplicates().values):

        df_group = df[df[group_column] == group]

        if not horizontal:
            axis.errorbar(x=df_group[f'{plot_entity}'], 
                        y=df_group['y'], 
                        xerr=df_group[[f'{plot_entity}_lower_distance', f'{plot_entity}_upper_distance']].values.T, 
                        label=group,
                        capsize=0, marker='s', markersize=markersize, 
                        markeredgecolor=colors[idx], 
                        markerfacecolor=colors[idx], linestyle='', elinewidth=1, color=colors[idx])
        else:
            axis.errorbar(y=df_group[f'{plot_entity}'], 
                        x=df_group['y'], 
                        yerr=df_group[[f'{plot_entity}_lower_distance', f'{plot_entity}_upper_distance']].values.T, 
                        label=group,
                        capsize=0, marker='s', markersize=markersize, 
                        markeredgecolor=colors[idx], 
                        markerfacecolor=colors[idx], linestyle='', elinewidth=1, color=colors[idx])


    # add labels on y-axis
    ticks = list(y_labels_dict.values())
    labels = list(y_labels_dict.keys())
    if not horizontal:
        axis.set_yticks(ticks)
        axis.set_yticklabels(labels)
    else:
        axis.set_xticks(ticks)
        axis.set_xticklabels(labels)

    # add legend
    if legend_loc:
        axis.legend(loc=legend_loc, numpoints=1, **legend_kwargs)

    return df


# def plot_errorbar_grouped_transposed(df, axis, y_column, group_column, 
#                                      legend_loc='lower right', legend_kwargs={}, group_scale=0.1, y_scale=1,
#                                      y_order=None, group_order=None, vertical_loc=1, colors=get_cmap("Accent").colors, 
#                                      plot_entity='odds_ratio'):
    
#     """
#         Makes a vertical errorbar plot for odds ratio, grouped by `group_column`.  
#     """
#     # we will use label encoder to calculate the positions
#     # of the data on the plot
#     le = preprocessing.LabelEncoder()

#     # calculate positions of the "groups"
#     if y_order is None:
#         df['y_group_loc'] = le.fit_transform(df[y_column])
#     else:
#         df['y_group_loc'] = df[y_column].apply(lambda x: y_order[::-1].index(x))

#     # calculate positions of each plot within "groups"
#     if group_order is None:
#         df['y'] = df['y_group_loc']*y_scale + le.fit_transform(df[group_column])*group_scale
#     else:
#         df['y'] = df['y_group_loc']*y_scale + df[group_column].apply(lambda x: group_order[::-1].index(x))*group_scale

#     # calculate the positions of labels on x-axis
#     y_labels_dict = df[[y_column, 'y']].groupby(y_column).mean().apply(lambda x: round(x, 2)).to_dict()['y']

#     # sort values in our df by position, so we control
#     # for the order in the legend
#     df = df.sort_values(by=['y'], ascending=True)

#     # plot each group errorbar separatly with its own color
#     for idx, group in enumerate(df[group_column].drop_duplicates().values):

#         df_group = df[df[group_column] == group]
        
#         axis.errorbar(y=df_group[f'{plot_entity}'], 
#                       x=df_group['y'], 
#                       yerr=df_group[[f'{plot_entity}_lower_distance', f'{plot_entity}_upper_distance']].values.T, 
#                       label=group,
#                       capsize=0, marker='s', markersize=3, markeredgecolor=colors[idx], 
#                       markerfacecolor=colors[idx], linestyle='', elinewidth=1, color=colors[idx])

#     # plot vertical line, denoting absence of effect 
#     axis.axhline(y=vertical_loc,  linestyle='--', color='salmon', linewidth=3)

#     # add labels on x-axis
#     axis.set_xticks(list(y_labels_dict.values()))
#     axis.set_xticklabels(list(y_labels_dict.keys()))

#     # add legend
#     if legend_loc:
#         axis.legend(loc=legend_loc, numpoints=1, **legend_kwargs)
    