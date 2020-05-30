import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import pickle
from scipy import stats



#################################################
"""
Source for the 2 Heaatmap functions
https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
"""
#################################################



def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):

    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    n_ticks = 20
    if data.shape[1]>n_ticks :
        ax.set_xticks(np.linspace(0 ,data.shape[1], num = n_ticks , dtype=int))
        ax.set_yticks(np.linspace(0 ,data.shape[0], num = n_ticks , dtype=int))
        ax.set_xticklabels(np.linspace(0 ,data.shape[1], num = n_ticks , dtype=int))
        ax.set_yticklabels(np.linspace(0 ,data.shape[0], num = n_ticks , dtype=int))
    else:
        ax.set_xticks(np.arange(data.shape[1]))
        ax.set_yticks(np.arange(data.shape[0]))
        # ... and label them with the respective list entries.
        ax.set_xticklabels(col_labels)
        ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-50, ha="right",rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=0.001)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    #https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html

    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


# LINEAR REGRESSION (LR)
# Plotting function for linear regression

def lr(a, b, y, sample, label_lr = "Normalized Linear regression ", label_sample = "Zipf random sample "):
    # a = starting point, b = end point, y = data,
    # if sample == 1 sampled data will be plotted with same slope as the one found in LR.
    # x values fo LR are defined in [a+1,b]

    x_LR = np.arange(1+a , b+1)
    # Why 1+a : y[0] is defined in x = 1, so if we choose a = 0 then x = 1.
    # Why b+1: np.arange exclude b, so if I need b included, it's necessary ro put b+1 inside np.arange
    y_LR = y[a:b]

    #scipy LR
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x_LR), np.log10(y_LR))

    print('LR m for real data =', slope)
    print('LR q for real data =', intercept)

    # LR values, eq.3 report
    c = 10 ** intercept
    y_LR = (np.arange(1, len(y_LR)+1) ** round(slope, 3)) * c

    #plotting linear regression
    ax1.plot(np.arange(1, len(y_LR)+1), y_LR/np.sum(y_LR), label= label_lr , linewidth=1.3, zorder=100, color ='orange')

    #plotting x_min Big Red Dot
    ax1.scatter(a, y[a], s=30, color='red', label = '$x_{min}=$'+str(a)+ ", $\\alpha$ = " + str(round(slope, 3)) )

    if sample == 1:
        # Sampled data
        random_zipf1 = stats.zipf.pmf(np.arange(1,b+1), a= -slope)
        # Sample sorted
        random_zipf1 =  - np.sort(-random_zipf1)
        random_zipf2 = random_zipf1[a:b+1]/np.sum(random_zipf1[a:b+1])
        print(np.sum(random_zipf1), len(random_zipf1), len(random_zipf2))
        # Plotting
        ax1.scatter(np.arange(1,len(random_zipf2)+1), random_zipf2, s=0.4, label=label_sample, color='green' )

    return y_LR[1:a+1]

#HEATMAPS plotting function.
def distance_heatmap(js_D, counties, title):
    itmap4 = plt.figure()
    ax4 = itmap4.add_subplot(1, 1, 1)
    im4, cbar4 = heatmap(js_D, counties, counties, ax=ax4, cmap='YlGn', cbarlabel=str(title))

    # Annotation of values on the squares, it is not good if the matrix is huge
    #texts = annotate_heatmap(im, valfmt="{x:.2f} ")

    fig.tight_layout()

#JS-KS RATIO plotting function
def ratio(ks,js):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(ks,js, s = 0.1)
    # Axis settings
    ax.set_xlabel('KS distance')
    ax.set_ylabel('JS distance')
    ax.grid(color='xkcd:light blue', linestyle='--', which='major', zorder=1, alpha=0.85)
    ax.grid(color='xkcd:light blue', linestyle='--', which='minor', zorder=1, alpha=0.35)
    ax.set_xticks(np.arange(0, 0.7, 0.05), minor=True)
    ax.set_yticks(np.arange(0, 0.4, 0.05), minor=True)
    #line
    x = np.arange(1,len(ks))
    #ax.plot(x,x, color = 'red')

#Which counties you want to be plotted? Also numpy array are ok.
#Chose a number between 1 and 3075
which_county = [1819] #[1, 12, 75, 1356, 3000]
xmin = 22

#One of the row files to extract some usefull informations
Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                    low_memory=False)

# total no. of counties
TOT_counties = len(Adata[1][:]) - 1

# List of total number of words per county, the order is the same as row data files
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                       index_col=[0])

# Empty array to be filled with Total number of words of the i-th county
word = np.zeros(3074)

# CALLING FIGURE
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

#IMPORTING DATA
with open(r'C:\Users\Salvo\Desktop\y.pkl', 'rb') as f:
    y = pickle.load(f)


#CALCULATIONG CDF
cdf_y = [[]]
for i in range (3076):
    cdf_y.append([])
for i in which_county:
    word[i] = int(tot_word.iloc[i, 1])
    #b is the total number of word for that county
    b = int(word[i])
    partial_sum = 0

    for j in range(b):
        partial_sum = partial_sum + y[i - 1][j]
        cdf_y[i].append(partial_sum)

    partial_sum = 0

    #PLOTTING
    print(str(i) + Adata[1][i]) # Plotting this county

    #PDF
    ax1.plot(np.arange(1,len(y[i-1][:b])+1),y[i-1][:b],
             linewidth = 1.3, label=str(i) + Adata[1][i] + str(b)+' words.', color='xkcd:navy blue')
    #CDF
    #ax1.plot(np.arange(1, len(y[i-1][:b]) + 1), cdf_y[i][:b],
    #         linewidth=0.5 ,label='CDF_'+str(i) + ' ' + Adata[1][i] + ', ' + str(b) + ' words.')

    #LINEAR REGRESSION plot
    y_LR = np.zeros(b+xmin)
    y_LR[1:xmin+1] = lr(xmin,b,y[i-1][:],1)
    y_LR[xmin:] = y[i - 1][:b]
    y_LR = y_LR/np.sum(y_LR)

    ax1.plot(np.arange(1, len(y_LR[1:]) + 1), y_LR[1:],
             linewidth=1.3, label='Modified PDF', color='brown')
    xxx = np.arange(1, 10000, 0.5)
    yyy = 1 / (5 + xxx) ** 1.5
    ax1.plot(xxx, yyy,
             linewidth=1.3, label='Modified Zipf', color='green')

# Figure settings

#x and y limits
ax1.set_xlim([0.1, 3*b ])
ax1.set_ylim([0.000000000001, 1])


#Axis grid
ax1.grid(color='xkcd:light blue', linestyle='--', which='major', zorder=1, alpha=0.85)
ax1.grid(color='xkcd:light blue', linestyle='--', which='minor', zorder=1, alpha=0.35)
ax1.set_xticks(np.arange(1,70000, 100), minor=True)
ax1.set_yticks(np.arange(1,70000, 100), minor=True)


# Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')

# Axis titles
ax1.set_xlabel('Ranking')
ax1.set_ylabel('Normalized frequency')

ax1.legend(loc='best')

########################################################################################################à
# IMPORTING DISTANCE MATRICES

#KS
#loading distance matrix
npzfile2 = np.load(r'C:\Users\Salvo\Desktop\Project\ks_distanceMatrix.npz')
print(npzfile2.files)
#assigning arrays in the file to new arrays
myD = npzfile2['ks_D']
#grid is the list of the counties in the order they appear in the ditsance matrix
counties1 = npzfile2['counties']

#JS
js_distance = np.load(r'C:\Users\Salvo\Desktop\Project\js_distanceMatrix.npz')
print(js_distance.files)
js_D = js_distance['js_D']
counties = js_distance['counties']

"""
#GEOGRAPHIC
combinations = 4723201
g_distance = np.load(r'C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\g_distances\\run\g_D_'+str(int(combinations))+'_combinations.npz')
print(g_distance.files)
g_D = g_distance['D']
g_grid = npzfile['counties']
print(g_grid)
"""

############
"""
PLOTTING HEATMAPS
"""
############

#KS Heatmap
#distance_heatmap(myD, counties1, 'KS heatmap')

#JS heatmap
#distance_heatmap(js_D, counties, 'JS heatmap')

#Geo heatmap
#distance_heatmap(g_D, g_grid, 'Geo heatmap')

#COMPARING DISTANCES, RATIO JS/KS PLOT
ratio(myD, js_D*(np.log2(10))**0.5)


plt.show()







