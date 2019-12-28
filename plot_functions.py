import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

def plot_range(xlabel, ylabel, title, x, values):
    """x and values should have the same size"""

    plt.plot(x, values, 'r-', linewidth=2)
    plt.gcf().set_size_inches(8,2);
    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel);

def plot_year_multi(*args):
    
    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(1, 1, 1)
    ox = np.arange(1,366,1)
    for arg in args:
        ax.plot(ox, arg[0], linewidth=2, label=arg[1])
        
    ax.set_ylabel(r'$values$'); ax.set_xlabel(r'$days$'); ax.legend(loc='best')
    
#def plot_year(xlabel, ylabel, title, values):
#    """values should be an array of 365 values"""
#
#    ox = np.arange(1,366,1)
#    plt.plot(ox, values, 'r-', linewidth=2)
#    plt.gcf().set_size_inches(8,2)
#    plt.title(title); plt.xlabel(xlabel); plt.ylabel(ylabel)

#def plot_triple(first, labelfirst, second, labelsecond, third, labelthird):
#    
#    # --- get an empty Figure and add an Axes
#    fig = plt.figure(figsize=(10, 4))
#    ax = fig.add_subplot(1, 1, 1) # row-col-num
#    # --- line plot data on the Axes
#    ox = np.arange(1,366,1)
#    ax.plot(ox, first, 'b-', linewidth=2, label=labelfirst)
#    ax.plot(ox, second, 'g-', linewidth=2, label=labelsecond)
#    ax.plot(ox, third, 'r-', linewidth=2, label=labelthird)
#    # --- add title and axis labels
#    #ax.set_title('The name')
#    ax.set_ylabel(r'$values$', fontsize=16)
#    ax.set_xlabel(r'$days$', fontsize=16)
#    # --- plot a legend in the best location
#    ax.legend(loc='best')
#    # --- add grid – not in default classic style
#    ax.grid(True)
#    # --- improve the layout
#    fig.tight_layout(pad=1)
    
if __name__ == '__main__':
    print('This is a plot functions module')
