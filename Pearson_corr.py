__author__ = "Vy Bui"
__email__ = "01bui@cua.edu"

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
import sys

def list_file_name(prefix, path):
    file_dir_list = []
    file_name_list = []
    for dirName, subdirList, fileList in os.walk(path):
        if subdirList == []:
            for name in fileList:
                if prefix in name:
                    file_dir_list.append(os.path.join(dirName, name))
                    file_name_list.append(name)
    return file_dir_list, file_name_list


def heatmap(data, ax=None,
            cbar_kw={}, cbarlabel="", grid_width=3, **kwargs):
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
    #ax.set_xticks(np.arange(data.shape[1]))
    #ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    #ax.set_xticklabels(col_labels)
    #ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    #ax.tick_params(top=True, bottom=False,
    #               labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    #plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
    #         rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=grid_width)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
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
            if data[i, j] >= 0.99:
                text = im.axes.text(j, i, 'N', **kw)
            else:
                text = im.axes.text(j, i, '', **kw)
            texts.append(text)

    return texts

def correlate(path1, path2):

    path = 'F:/CellPaintingAssay/FeatureMeasurement_CellProfiler/U2OS_CPA_FM_2_and_3_Predict_Final'
    MeasureCorrelation_pre, _ = list_file_name('MeasureCorrelation', path1)
    MeasureGranularity_pre, _ = list_file_name('MeasureGranularity', path1)
    MeasureObjectIntensity_pre, _ = list_file_name('MeasureObjectIntensity', path1)
    MeasureObjectRadialDistribution_pre, _ = list_file_name('MeasureObjectRadialDistribution', path1)
    MeasureObjectSizeShape_pre, _ = list_file_name('MeasureObjectSizeShape', path1)
    MeasureTexture_pre, _ = list_file_name('MeasureTexture', path1)

    MeasureCorrelation_GT, _ = list_file_name('MeasureCorrelation', path2)
    MeasureGranularity_GT, _ = list_file_name('MeasureGranularity', path2)
    MeasureObjectIntensity_GT, _ = list_file_name('MeasureObjectIntensity', path2)
    MeasureObjectRadialDistribution_GT, _ = list_file_name('MeasureObjectRadialDistribution', path2)
    MeasureObjectSizeShape_GT, _ = list_file_name('MeasureObjectSizeShape', path2)
    MeasureTexture_GT, _ = list_file_name('MeasureTexture', path2)


    # MeasureCorrelation ============================================================================
    print("MeasureCorrelation - calculation")
    ROW = 120
    MeasureCorrelation_pre_array = np.zeros((ROW, len(MeasureCorrelation_pre)), dtype=np.float)
    MeasureCorrelation_GT_array = np.zeros((ROW, len(MeasureCorrelation_GT)), dtype=np.float)
    for idx_img in range(len(MeasureCorrelation_GT)):
        MeasureCorrelation_pre_value = pd.read_csv(MeasureCorrelation_pre[idx_img])
        MeasureCorrelation_pre_array[:, idx_img] = MeasureCorrelation_pre_value['Value'].values

        MeasureCorrelation_GT_value = pd.read_csv(MeasureCorrelation_GT[idx_img])
        MeasureCorrelation_GT_array[:, idx_img] = MeasureCorrelation_GT_value['Value'].values

    R_corr_MeasureCorrelation = np.zeros(ROW, dtype=np.float)
    for i in range(ROW):
        R_corr_MeasureCorrelation[i] = np.corrcoef(MeasureCorrelation_pre_array[i], MeasureCorrelation_GT_array[i])[0, 1]

    # MeasureGranularity ============================================================================
    print("MeasureGranularity - calculation")
    ROW = 5
    COL = 16
    MeasureGranularity_pre_array = np.zeros((ROW, COL, len(MeasureGranularity_pre)), dtype=np.float)
    MeasureGranularity_GT_array = np.zeros((ROW, COL, len(MeasureGranularity_pre)), dtype=np.float)

    for idx_img in range(len(MeasureGranularity_GT)):

        GS_pre = np.zeros((ROW, COL), dtype=np.float)
        GS_GT = np.zeros((ROW, COL), dtype=np.float)
        for col in range(COL):
            name = 'GS' +str(col + 1)

            MeasureGranularity_pre_value = pd.read_csv(MeasureGranularity_pre[idx_img])
            GS_pre[:, col] = MeasureGranularity_pre_value[name].values

            MeasureGranularity_GT_value = pd.read_csv(MeasureGranularity_GT[idx_img])
            GS_GT[:, col] = MeasureGranularity_GT_value[name].values

        MeasureGranularity_pre_array[:, :, idx_img] = GS_pre
        MeasureGranularity_GT_array[:, :, idx_img] = GS_GT

    R_corr_MeasureGranularity = np.zeros((ROW, 16), dtype=np.float)
    for rol in range(ROW):
        for col in range(COL):
            R_corr_MeasureGranularity[rol, col] = np.corrcoef(MeasureGranularity_pre_array[rol, col, :],
                                              MeasureGranularity_GT_array[rol, col, :])[0, 1]

    # MeasureObjectIntensity ============================================================================
    print("MeasureObjectIntensity - calculation")
    ROW = 285
    COL = 3
    MeasureObjectIntensity_pre_array = np.zeros((ROW, COL, len(MeasureObjectIntensity_pre)), dtype=np.float)
    MeasureObjectIntensity_GT_array = np.zeros((ROW, COL, len(MeasureObjectIntensity_GT)), dtype=np.float)

    for idx_img in range(len(MeasureObjectIntensity_GT)):
        MeasureObjectIntensity_pre_value = pd.read_csv(MeasureObjectIntensity_pre[idx_img])
        MeasureObjectIntensity_pre_array[:, 0, idx_img] = MeasureObjectIntensity_pre_value['Mean'].values
        MeasureObjectIntensity_pre_array[:, 1, idx_img] = MeasureObjectIntensity_pre_value['Median'].values
        MeasureObjectIntensity_pre_array[:, 2, idx_img] = MeasureObjectIntensity_pre_value['STD'].values

        MeasureObjectIntensity_GT_value = pd.read_csv(MeasureObjectIntensity_GT[idx_img])
        MeasureObjectIntensity_GT_array[:, 0, idx_img] = MeasureObjectIntensity_GT_value['Mean'].values
        MeasureObjectIntensity_GT_array[:, 1, idx_img] = MeasureObjectIntensity_GT_value['Median'].values
        MeasureObjectIntensity_GT_array[:, 2, idx_img] = MeasureObjectIntensity_GT_value['STD'].values

    R_corr_MeasureObjectIntensity = np.zeros((ROW, 3), dtype=np.float)
    for row in range(ROW):
        for col in range(COL):
            R_corr_MeasureObjectIntensity[row, col] = np.corrcoef(MeasureObjectIntensity_pre_array[row, col, :],
                                                                  MeasureObjectIntensity_GT_array[row, col, :])[0, 1]

    # MeasureObjectRadialDistribution ============================================================================
    print("MeasureObjectRadialDistribution - calculation")
    ROW = 48
    COL = 3
    MeasureObjectRadialDistribution_pre_array = np.zeros((ROW, COL, len(MeasureObjectRadialDistribution_pre)), dtype=np.float)
    MeasureObjectRadialDistribution_GT_array = np.zeros((ROW, COL, len(MeasureObjectRadialDistribution_GT)), dtype=np.float)

    for idx_img in range(len(MeasureObjectRadialDistribution_GT)):
        MeasureObjectRadialDistribution_pre_value = pd.read_csv(MeasureObjectRadialDistribution_pre[idx_img])
        MeasureObjectRadialDistribution_pre_array[:, 0, idx_img] = MeasureObjectRadialDistribution_pre_value['Fraction'].values
        MeasureObjectRadialDistribution_pre_array[:, 1, idx_img] = MeasureObjectRadialDistribution_pre_value['Intensity'].values
        MeasureObjectRadialDistribution_pre_array[:, 2, idx_img] = MeasureObjectRadialDistribution_pre_value['COV'].values

        MeasureObjectRadialDistribution_GT_value = pd.read_csv(MeasureObjectRadialDistribution_GT[idx_img])
        MeasureObjectRadialDistribution_GT_array[:, 0, idx_img] = MeasureObjectRadialDistribution_GT_value['Fraction'].values
        MeasureObjectRadialDistribution_GT_array[:, 1, idx_img] = MeasureObjectRadialDistribution_GT_value['Intensity'].values
        MeasureObjectRadialDistribution_GT_array[:, 2, idx_img] = MeasureObjectRadialDistribution_GT_value['COV'].values

    R_corr_MeasureObjectRadialDistribution = np.zeros((ROW, 3), dtype=np.float)
    for row in range(ROW):
        for col in range(COL):
            R_corr_MeasureObjectRadialDistribution[row, col] = np.corrcoef(MeasureObjectRadialDistribution_pre_array[row, col, :],
                                                                  MeasureObjectRadialDistribution_GT_array[row, col, :])[0, 1]

    # MeasureObjectSizeShape ============================================================================
    print("MeasureObjectSizeShape - calculation")
    ROW = 54  # No Zernike
    COL = 3

    MeasureObjectSizeShape_pre_value = pd.read_csv(MeasureObjectSizeShape_pre[0])
    Feature = MeasureObjectSizeShape_pre_value['Feature'].values

    MeasureObjectSizeShape_pre_array = np.zeros((ROW, COL, len(MeasureObjectSizeShape_pre)), dtype=np.float)
    MeasureObjectSizeShape_GT_array = np.zeros((ROW, COL, len(MeasureObjectSizeShape_GT)), dtype=np.float)

    for idx_img in range(len(MeasureObjectSizeShape_GT)):
        MeasureObjectSizeShape_pre_value = pd.read_csv(MeasureObjectSizeShape_pre[idx_img])
        MeasureObjectSizeShape_GT_value = pd.read_csv(MeasureObjectSizeShape_GT[idx_img])

        MeasureObjectSizeShape_pre_NoZernike_Mean = []
        MeasureObjectSizeShape_pre_NoZernike_Median = []
        MeasureObjectSizeShape_pre_NoZernike_STD = []

        MeasureObjectSizeShape_GT_NoZernike_Mean = []
        MeasureObjectSizeShape_GT_NoZernike_Median = []
        MeasureObjectSizeShape_GT_NoZernike_STD = []

        for i in range(len(Feature)):
            if 'Zernike' not in Feature[i]:
                MeasureObjectSizeShape_pre_NoZernike_Mean.append(MeasureObjectSizeShape_pre_value['Mean'].values[i])
                MeasureObjectSizeShape_pre_NoZernike_Median.append(MeasureObjectSizeShape_pre_value['Median'].values[i])
                MeasureObjectSizeShape_pre_NoZernike_STD.append(MeasureObjectSizeShape_pre_value['Std'].values[i])

                MeasureObjectSizeShape_GT_NoZernike_Mean.append(MeasureObjectSizeShape_GT_value['Mean'].values[i])
                MeasureObjectSizeShape_GT_NoZernike_Median.append(MeasureObjectSizeShape_GT_value['Median'].values[i])
                MeasureObjectSizeShape_GT_NoZernike_STD.append(MeasureObjectSizeShape_GT_value['Std'].values[i])

        MeasureObjectSizeShape_pre_array[:, 0, idx_img] = MeasureObjectSizeShape_pre_NoZernike_Mean
        MeasureObjectSizeShape_pre_array[:, 1, idx_img] = MeasureObjectSizeShape_pre_NoZernike_Median
        MeasureObjectSizeShape_pre_array[:, 2, idx_img] = MeasureObjectSizeShape_pre_NoZernike_STD

        MeasureObjectSizeShape_GT_array[:, 0, idx_img] = MeasureObjectSizeShape_GT_NoZernike_Mean
        MeasureObjectSizeShape_GT_array[:, 1, idx_img] = MeasureObjectSizeShape_GT_NoZernike_Median
        MeasureObjectSizeShape_GT_array[:, 2, idx_img] = MeasureObjectSizeShape_GT_NoZernike_STD

    R_corr_MeasureObjectSizeShape = np.zeros((ROW, COL), dtype=np.float)
    for row in range(ROW):
        for col in range(COL):
            R_corr_MeasureObjectSizeShape[row, col] = np.corrcoef(MeasureObjectSizeShape_pre_array[row, col, :],
                                                                  MeasureObjectSizeShape_GT_array[row, col, :])[0, 1]

    # MeasureTexture ============================================================================
    print("MeasureTexture - calculation")
    ROW = 3360
    MeasureTexture_pre_array = np.zeros((ROW, len(MeasureTexture_pre)), dtype=np.float)
    MeasureTexture_GT_array = np.zeros((ROW, len(MeasureTexture_GT)), dtype=np.float)
    for idx_img in range(len(MeasureTexture_GT)):
        MeasureTexture_pre_value = pd.read_csv(MeasureTexture_pre[idx_img])
        MeasureTexture_pre_array[:, idx_img] = MeasureTexture_pre_value['Value'].values

        MeasureTexture_GT_value = pd.read_csv(MeasureTexture_GT[idx_img])
        MeasureTexture_GT_array[:, idx_img] = MeasureTexture_GT_value['Value'].values

    R_corr_MeasureTexture = np.zeros(ROW, dtype=np.float)
    for i in range(ROW):
        R_corr_MeasureTexture[i] = np.corrcoef(MeasureTexture_pre_array[i], MeasureTexture_GT_array[i])[0, 1]


    A = np.reshape(R_corr_MeasureCorrelation, (24, 5))
    B = np.reshape(R_corr_MeasureObjectIntensity, (57, 15))
    C = np.reshape(R_corr_MeasureObjectRadialDistribution, (24, 6))
    D = np.reshape(R_corr_MeasureObjectSizeShape, (18, 9))
    E = np.reshape(R_corr_MeasureTexture, (120, 28))

    return A, R_corr_MeasureGranularity, B, C, D, E


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise IOError("Not enough inputs")

    path1 = sys.argv[1]
    path2 = sys.argv[2]

    A, B, C, D, E, F = correlate(path1, path2)

    fig, ax = plt.subplots()
    im, cbar = heatmap(A, ax=ax, cmap="magma_r", grid_width=3)
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    im, cbar = heatmap(B, ax=ax, cmap="magma_r", grid_width=3)
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    im, cbar = heatmap(C, ax=ax, cmap="magma_r", grid_width=1)
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    im, cbar = heatmap(D, ax=ax, cmap="magma_r", grid_width=3)
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    im, cbar = heatmap(E, ax=ax, cmap="magma_r")
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots()
    im, cbar = heatmap(F, ax=ax, cmap="magma_r", grid_width=1)
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.show()



