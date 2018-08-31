import numpy as np
import matplotlib.pyplot as plt

def plot_over_year(xlabel, ylabel, title, values):
    """values should be an array of 365 values"""

    # --- Select a style
    plt.style.use('classic')
    
    ox = np.arange(1,366,1)
    plt.plot(ox, values);
    plt.gcf().set_size_inches(10,2);
    plt.title(title);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);

def plot_over_range(xlabel, ylabel, title, ox, values):
    """ox and values should have the same size"""
    
    plt.plot(ox, values);
    plt.gcf().set_size_inches(10,2);
    plt.title(title);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);

def plot_example():
    # --- get an empty Figure and add an Axes
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(1, 1, 1) # row-col-num
    # --- line plot data on the Axes
    ox = np.arange(1,366,1)
    ax.plot(ox, no3_limiter, 'b-', linewidth=2, label=r'$NO_{3}^{-}$ limiter')
    ax.plot(ox, po4_limiter, 'g-', linewidth=2, label=r'$PO_{4}^{2-}$ limiter')
    ax.plot(ox, si_limiter, 'r-', linewidth=2, label=r'$Si$ limiter')
    # --- add title and axis labels
    ax.set_title('The limiters')
    ax.set_ylabel(r'$limiters$', fontsize=16)
    ax.set_xlabel(r'$days$', fontsize=16)
    # --- plot a legend in the best location
    ax.legend(loc='best')
    # --- add grid – not in default classic style
    ax.grid(True)
    # --- improve the layout
    fig.tight_layout(pad=1)

if __name__ == '__main__':
    print('This is a plot functions module')