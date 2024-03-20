import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import matplotlib 
import numpy as np
import shap

def reorderLegend(ax=None,order=None,unique=False):
    if ax is None: ax=plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0])) # sort both labels and handles by labels
    if order is not None: # Sort according to a given list (not necessarily complete)
        keys=dict(zip(order,range(len(order))))
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t,keys=keys: keys.get(t[0],np.inf)))
    if unique:  labels, handles= zip(*unique_everseen(zip(labels,handles), key = labels)) # Keep only the first of each handle
    ax.legend(handles, labels, loc = 'lower right')
    return(handles, labels)

def unique_everseen(seq, key=None):
    seen = set()
    seen_add = seen.add
    return [x for x,k in zip(seq,key) if not (k in seen or seen_add(k))]

def plot_prediction(y_test, list_of_predictions, label, save_prefix = ''):
    matplotlib.style.use('seaborn')

    # n = len(sorted_values)

    fig, axs = plt.subplots(3, 3, figsize=(9, 9))
    for i, ax in enumerate(axs.flatten()):
        for j in range(1):
            index = i * 1 + j
            sorted_values = [pred[y_test[index].reset_index()['Age_death'].sort_values().index] for pred in list_of_predictions]
            ax.plot(sorted_values[index], label=f'{label[index]}')
            ax.legend(loc = 'upper left' )
        ax.plot(np.array(y_test[index].reset_index()['Age_death'].sort_values()), label = 'Ground Truth')
        ax.legend(loc = 'upper left' )
        ax.set_xlabel('Index')
        ax.set_ylabel('Predicted Age')
    
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

    matplotlib.rc('font', **font)
    plt.tight_layout()
    time = datetime.datetime.now().strftime('%Y-%m-%d')
    plt.savefig(f'YL_{time}_{save_prefix}_prediction.pdf')
    plt.show()
    plt.close()


def plot_shap(prefix, shap_values, sort_by_age = True):
    # size of the figure
    plt.figure(figsize=(15, 10))
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    shap.plots.heatmap(shap_values, shap_values.sum(1), show=False) if sort_by_age else shap.plots.heatmap(shap_values, show=False)
    plt.savefig(f'{prefix}_shap.pdf')

def plot_table(ret, save_prefix = ''):
    font = {'size'   : 15}
    matplotlib.rc('font', **font)

    ret = ret.set_index('Celltype')
    ret = ret.sort_values('Spearman Correlation')
    fig, ax = plt.subplots(figsize=(15, 9))
    ret[['Pearson Correlation', 'Spearman Correlation', ]].plot(kind='barh', color=['skyblue', 'coral'], ax=ax)
    ax.legend()
    ax.set_xlabel('Correlation')
    sample_num = list(set(ret['number of samples']))
    gene_num = list(set(ret['number of genes']))
    sample_num_str = f'Fixed to {sample_num[0]}' if len(sample_num) == 1 else 'Unfixed' 
    ax.set_title(f'Spearman and Pearson Correlation, DatasetID {save_prefix}, Sample Size {sample_num_str}, Gene number {gene_num[0]}')
    ax.grid(axis='x')
    # Add numeric labels to the bars
    for container in ax.containers:
        ax.bar_label(container, fmt='%.2f')

    # Add the number of samples as text annotations next to the cell types
    if len(sample_num) != 1:
        for i, cell_type in enumerate(ret.index):
            num_samples = ret.loc[cell_type, "number of samples"]
            ax.text(1, i, f' ({num_samples})', va='center', ha='left', transform=ax.get_yaxis_transform(), color='red')

    reorderLegend(ax,['Spearman Correlation', 'Pearson Correlation'])
    # Save the plot as a PDF file
    # get today's date
    date = datetime.datetime.now().strftime('%Y-%m-%d')
    str_suffix = f'Sample_Size_{sample_num_str}_Gene_number_{gene_num[0]}'
    plt.savefig(f'YL_{save_prefix}_{date}_{str_suffix}.pdf', format='pdf')
    plt.close()
# plt.close()