import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import pickle


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
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",rotation_mode="anchor")

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
#One of the row files to extract some usefull informations
Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                    low_memory=False)

# total no. of counties
TOT_counties = len(Adata[1][:]) - 1

# List of total number of words per county, the order is the same as row data files
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                       index_col=[0])
print(tot_word)

#Which counties you want to be plotted? Also numpy array are ok.
#Chose a number between 1 and 3075
which_county = [1, 141, 12, 200, 134, 103, 34, 84, 100, 180]


# Empty array to be filled with Total number of words of the i-th county
word = np.zeros(3074)

#c is the index to identify the i-th county
c = which_county

#I used this part when I needed to produce every time the list for y and cdf_y.
# Now I have a file for y and better code for cdf_y
"""
# Initializing variables
x = np.arange(1, 70000)
for i in  range (len(which_county)) :
    word[i] = int(tot_word.iloc[c[i], 1])
    b = int(word[i])

    y = [[0 for indiceacaso in range(1)] for indiceacaso2 in which_county]
    cdf_y = [[0 for indiceacaso in range(70000)] for indiceacaso2 in which_county]
"""

# GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

#IMPORTING DATA
with open(r'C:\Users\Salvo\Desktop\y.pkl', 'rb') as f:
    y = pickle.load(f)

cdf_y = [[]]
#CDF
for i in which_county:
    word[i] = int(tot_word.iloc[c[i], 1])
    b = int(word[i])
    cdf_y.append([])
    partial_sum = 0
    # print(b, len(y[i]))

    for j in range(b):
        partial_sum = partial_sum + y[i - 1][j]
        cdf_y[i].append(partial_sum)


    partial_sum = 0
    print(cdf_y[i])


    #GRAPHIC y and cdf

    ax1.plot(np.arange(1,len(y[i][:b])+1),y[i][:b],
             linewidth = 0.5, label='PMF_'+str(c[i])+' '+ Adata[1][c[i]] +', '+str(b)+' words.')

    ax1.plot(np.arange(1, len(y[i][:b]) + 1), cdf_y[i][:b],
             linewidth=0.5 ,label='CDF_'+str(c[i]) + ' ' + Adata[1][c[i]] + ', ' + str(b) + ' words.')


#Axis grid
ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
ax1.set_xticks(np.arange(1,10000, 100), minor=False)



# Figure settings

# Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')

# Axis titles
ax1.set_xlabel('Ranking')
ax1.set_ylabel('Normalized frequency')

ax1.legend(loc='best')



############
"""
HEATMAP
"""
##########



#######  KS Heatmaps

#calling fig object
itmap = plt.figure()
ax = itmap.add_subplot(1, 1, 1)

#loading distance matrix
combinations = 19900
npzfile = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\ks\run\myD_value_'+str(int(combinations))+'_combinations.npz')
print(npzfile.files)

#assigning arrays in the file to new arrays
D = npzfile['D']
#p = npzfile['p']

#grid is the list of the counties in the order they appear in the ditsance matrix
grid = npzfile['counties']
print(grid)

#the same as before
combinations = 4723201
npzfile2 = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\ks\run\myD_value_'+str(int(combinations))+'_combinations.npz')
print(npzfile2.files)
myD = npzfile2['D']
grid2 = npzfile2['counties']



#calling plot functions
im, cbar = heatmap(D, grid, grid, ax=ax, cmap="YlGn", cbarlabel="KS distance "+str(grid[0])+'-'+str(grid[len(grid)-1]))
itmap3 = plt.figure()

ax3 = itmap3.add_subplot(1, 1, 1)
im3, cbar3 = heatmap(myD, grid, grid, ax=ax3, cmap="YlGn", cbarlabel="myKS distance "+str(grid[0])+'-'+str(grid[len(grid)-1]))


#Annotation of values on the squares, it is not good if the matrix is huge
#texts = annotate_heatmap(im, valfmt="{x:.2f} ")
#fig.tight_layout()


"""
#G Heatmap
itmap2 = plt.figure()
ax2 = itmap2.add_subplot(1, 1, 1)

combinations = 4723201
g_distance = np.load(r'C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\g_distances\\run\g_D_'+str(int(combinations))+'_combinations.npz')
print(g_distance.files)
g_D = g_distance['D']
g_grid = npzfile['counties']
print(g_grid)
im, cbar = heatmap(g_D, grid, grid, ax=ax2, cmap="YlGn", cbarlabel="Geographic distance"+str(grid[0])+'-'+str(grid[len(grid)-1]))
#texts = annotate_heatmap(im, valfmt="{x:.2f} ")
fig.tight_layout()
"""


#JS heatmap
itmap4 = plt.figure()
ax4 = itmap4.add_subplot(1, 1, 1)

combinations = 4723201
js_distance = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\js\run\jsD_value_'+str(int(combinations))+'_combinations.npz')
print(js_distance.files)
js_D = js_distance['jd_D']
js_grid = npzfile['counties']
#print(js_grid)

im4, cbar4 = heatmap(js_D, grid, grid, ax=ax4, cmap="YlGn", cbarlabel="JS distance"+str(grid[0])+'-'+str(grid[len(grid)-1]))
#texts = annotate_heatmap(im, valfmt="{x:.2f} ")
fig.tight_layout()


#COMPARING DISTANCES
figDD = plt.figure()
axD = figDD.add_subplot(1, 1, 1)
axD.scatter(myD, js_D, s = 0.1)
# Axis titles
axD.set_xlabel('myKS distance')
axD.set_ylabel('js distance')

plt.show()







