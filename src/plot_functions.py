import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def plot_range(xlabel, ylabel, title, x, values):
    """x and values should have the same size"""

    plt.plot(x, values, 'r-', linewidth=2)
    plt.gcf().set_size_inches(8, 2)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


def plot_year_multi(*args):
    """*args should be iterateble with 2 elements (value, label);
       value should be an array of 365 elements"""

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(1, 1, 1)
    ox = np.arange(1, 366, 1)
    for arg in args:
        ax.plot(ox, arg[0], linewidth=2, label=arg[1])

    ax.set_ylabel(r'$values$')
    ax.set_xlabel(r'$days$')
    ax.legend(loc='best')


if __name__ == '__main__':
    print('This is a plot functions module')
