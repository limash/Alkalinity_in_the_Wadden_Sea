# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 0.8.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import numpy as np
import matplotlib.pyplot as plt

def plot_over_year(xlabel, ylabel, title, values):
    """values should be an array of 365 values"""

    # --- Select a style
    plt.style.use('classic')

    ox = np.arange(1,366,1)
    plt.plot(ox, values);
    plt.gcf().set_size_inches(12,4);
    plt.title(title);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);

def plot_over_range(xlabel, ylabel, title, ox, values):
    """ox and values should have the same size"""

    # --- Select a style
    plt.style.use('classic')

    plt.plot(ox, values);
    plt.gcf().set_size_inches(12,2);
    plt.title(title);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);

def plot_triple(first, labelfirst, second, labelsecond, third, labelthird):

    # --- Select a style
    plt.style.use('classic')
    
    # --- get an empty Figure and add an Axes
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    # --- line plot data on the Axes
    ox = np.arange(1,366,1)
    ax.plot(ox, first, 'b-', linewidth=2, label=labelfirst)
    ax.plot(ox, second, 'g-', linewidth=2, label=labelsecond)
    ax.plot(ox, third, 'r-', linewidth=2, label=labelthird)
    # --- add title and axis labels
    #ax.set_title('The name')
    ax.set_ylabel(r'$values$', fontsize=16)
    ax.set_xlabel(r'$days$', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='best')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)

def plot_double(first, labelfirst, second, labelsecond):

    # --- Select a style
    plt.style.use('classic')
    
    # --- get an empty Figure and add an Axes
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    # --- line plot data on the Axes
    ox = np.arange(1,366,1)
    ax.plot(ox, first, 'b-', linewidth=2, label=labelfirst)
    ax.scatter(ox, second, c='r', label=labelsecond)
    # --- add title and axis labels
    #ax.set_title('The name')
    ax.set_ylabel(r'$values$', fontsize=16)
    ax.set_xlabel(r'$days$', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='best')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)
if __name__ == '__main__':
    print('This is a plot functions module')
