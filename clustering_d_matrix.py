import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt

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
    ax.grid(which="minor", color="w", linestyle='-', linewidth=0.0005)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar



#KS
combinations = 4723201
myks_file = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\ks\run\myD_value_'+str(int(combinations))+'_combinations.npz')
myD = myks_file['D']
ks_grid = myks_file['counties']


#JS
combinations = 4723201
js_distance = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\js\run\jsD_value_'+str(int(combinations))+'_combinations.npz')
print(js_distance.files)
js_D = js_distance['jd_D']
js_grid = js_distance['counties']





# generate the linkage matrix
Z = linkage(myD)
Z_2 = linkage(js_D)
print(myD)

# calculate full dendrogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram, KS 1')
plt.xlabel('sample index')
plt.ylabel('distance')

dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=6.,  # font size for the x axis labels
    truncate_mode='lastp',
    p=180,
    #show_contracted=True,
      # useful in small plots so annotations don't overlap
)



# calculate full dendrogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram, KS 2')
plt.xlabel('sample index')
plt.ylabel('distance')

dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=6.,  # font size for the x axis labels
    truncate_mode='lastp',
    p=270,
    #show_contracted=True,
      # useful in small plots so annotations don't overlap
)



# calculate full dendrogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram, JS 1')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    Z_2,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=6.,  # font size for the x axis labels
    truncate_mode='lastp',
    p=360,
    show_contracted=True,
    show_leaf_counts = True,
    count_sort = False,


      # useful in small plots so annotations don't overlap
)

plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram, JS 2')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    Z_2,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=6.,  # font size for the x axis labels
    truncate_mode='lastp',
    p=250,
    #show_contracted=True,
      # useful in small plots so annotations don't overlap
)


itmap = plt.figure()
ax = itmap.add_subplot(1, 1, 1)
im, cbar = heatmap(myD, ks_grid, ks_grid, ax=ax, cmap="YlGn", cbarlabel="KS distance "+str(ks_grid[0])+'-'+str(ks_grid[len(ks_grid)-1]))

itmap4 = plt.figure()
ax4 = itmap4.add_subplot(1, 1, 1)
im4, cbar4 = heatmap(js_D, js_grid, js_grid, ax=ax4, cmap="YlGn", cbarlabel="JS distance"+str(js_grid[0])+'-'+str(js_grid[len(js_grid)-1]))


bigV_ks = []
bigV_js = []

'''
for i in range (0, len(ks_grid)):
    for j in range (i+1, len(ks_grid)):
        if myD[i][j]>0.3:
            a = (i, j, myD[i][j])
            bigV_ks.append(a)
        if js_D[i][j]>0.3:
            a = (i, j, js_D[i][j])
            bigV_js.append(a)

print('ks: ',bigV_ks)
print('js: ',bigV_js)
'''


#Knee plot
fug2 = plt.figure()
ax2 = fug2.add_subplot(1, 1, 1)




#KS Knee method
last = Z[:, 2]
#print('last: ', last)
last_revERSED = last[::-1]
#print(last_revERSED)
idxs = np.arange(1, len(last) + 1)
ax2.plot(idxs, last_revERSED, linewidth = 0.5, label = 'ks distance values')

acceleration = np.diff(last, 2)  # 2nd derivative of the distances
acceleration_rev = acceleration[::-1]
ax2.plot(idxs[:-2] + 1, acceleration_rev, linewidth = 0.5, label = 'ks d acceleration')

#Number of cluster given from max of acceleration
k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters
print ("clusters:", k)




#JS Knee method
last = Z_2[:, 2]
print('last: ', last)
last_revERSED = last[::-1]
print(last_revERSED)
idxs = np.arange(1, len(last) + 1)
ax2.plot(idxs, last_revERSED, linewidth = 0.5, label = 'js distance values')

#Acceleration
acceleration = np.diff(last, 2)  # 2nd derivative of the distances
acceleration_rev = acceleration[::-1]
ax2.plot(idxs[:-2] + 1, acceleration_rev, linewidth = 0.5, label = 'js d acceleration')

#Number of cluster given from max of acceleration
k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters
print ("clusters:", k)


plt.xlabel('iteration')

plt.show()
