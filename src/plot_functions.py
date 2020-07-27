import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xr
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


def extract_alk(data_train):
    ds = xr.open_dataset(data_train[0])
    alk_df = ds['B_C_Alk'].to_dataframe()
    alk_surface = alk_df.groupby('z').get_group(data_train[1])
    alk = alk_surface.loc['2011-01-01':'2011-12-31']
    alk = alk.reset_index()
    return alk


def show_alk(data_train):
    fig = plt.figure(figsize=(10, 2))
    ax = fig.add_subplot(1, 1, 1)
    for item in data_train:
        ax.plot(item[0]['time'],
                item[0]['B_C_Alk'], linewidth=2, label=item[1])
    ax.legend(loc='best')
    ax.set_title('Alkalinity in the surface layer')
    plt.show()
 

if __name__ == '__main__':
    print('This is a plot functions module')