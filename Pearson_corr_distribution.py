__author__ = "Vy Bui"
__email__ = "01bui@cua.edu"

import os
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from Pearson_corr import heatmap, annotate_heatmap
import matplotlib.pyplot as plt

def parse_image_measurement(path):
    list_dir = os.listdir(path)
    all_cells = []
    all_nuclei = []
    all_cytoplasm = []
    for dir in list_dir:
        list_file = os.listdir(path + '/' + dir)

        Cells = []
        Nuclei = []
        Cytoplasm = []

        for file in list_file:
            if file.endswith('Image_MeasureGranularity.csv'):
                continue
            elif file.endswith('Image_MeasureCorrelation.csv'):
                #print(path + '/' + dir + '/' + file)
                sheet = pd.read_csv(path + '/' + dir+'/'+file, dtype=str)
                data = sheet.loc[:, ['Object', 'Value']]
                data = data.set_index("Object", drop=False)
                data.Value = data.Value.astype(float)

                if data.index.unique().__contains__('Cells'):
                    Cells.append(np.asarray(data.loc["Cells", "Value"], dtype=np.float))
                elif data.index.unique().__contains__('Cell'):
                    Cells.append(np.asarray(data.loc["Cell", "Value"], dtype=np.float))

                Nuclei.append(np.asarray(data.loc["Nuclei", "Value"], dtype=np.float))

                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Value"], dtype=np.float))
            elif file.endswith('Image_MeasureObjectIntensity.csv'):
                #print(path+ '/' + dir+'/'+file)
                sheet = pd.read_csv(path + '/' + dir + '/' + file, dtype=str)
                data = sheet.loc[:, ['Objects', 'Mean', 'Median', 'STD']]
                data = data.set_index("Objects", drop=False)
                data.Mean = data.Mean.astype(float)
                data.Median = data.Median.astype(float)
                data.STD = data.STD.astype(float)

                if data.index.unique().__contains__('Cells'):
                    Cells.append(np.asarray(data.loc["Cells", "Mean"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cells", "Median"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cells", "STD"], dtype=np.float))
                elif data.index.unique().__contains__('Cell'):
                    Cells.append(np.asarray(data.loc["Cell", "Mean"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cell", "Median"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cell", "STD"], dtype=np.float))

                Nuclei.append(np.asarray(data.loc["Nuclei", "Mean"], dtype=np.float))
                Nuclei.append(np.asarray(data.loc["Nuclei", "Median"], dtype=np.float))
                Nuclei.append(np.asarray(data.loc["Nuclei", "STD"], dtype=np.float))

                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Mean"], dtype=np.float))
                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Median"], dtype=np.float))
                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "STD"], dtype=np.float))
            elif file.endswith('Image_MeasureObjectRadialDistribution.csv'):
                #print(path+'/' + dir+'/'+file)
                sheet = pd.read_csv(path + '/' + dir + '/' + file, dtype=str)
                data = sheet.loc[:, ['Object', 'Fraction', 'Intensity', 'COV']]
                data = data.set_index("Object", drop=False)
                data.Fraction = data.Fraction.astype(float)
                data.Intensity = data.Intensity.astype(float)
                data.COV = data.COV.astype(float)

                if data.index.unique().__contains__('Cells'):
                    Cells.append(np.asarray(data.loc["Cells", "Fraction"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cells", "Intensity"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cells", "COV"], dtype=np.float))
                elif data.index.unique().__contains__('Cell'):
                    Cells.append(np.asarray(data.loc["Cell", "Fraction"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cell", "Intensity"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cell", "COV"], dtype=np.float))

                Nuclei.append(np.asarray(data.loc["Nuclei", "Fraction"], dtype=np.float))
                Nuclei.append(np.asarray(data.loc["Nuclei", "Intensity"], dtype=np.float))
                Nuclei.append(np.asarray(data.loc["Nuclei", "COV"], dtype=np.float))

                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Fraction"], dtype=np.float))
                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Intensity"], dtype=np.float))
                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "COV"], dtype=np.float))

            elif file.endswith('Image_MeasureObjectSizeShape.csv'):
                #print(path+'/' + dir+'/'+file)
                sheet = pd.read_csv(path + '/' + dir + '/' + file, dtype=str)
                data = sheet.loc[:, ['Objects', 'Feature', 'Mean', 'Median', 'Std']]
                data = data.set_index("Feature", drop=False)
                data = data.loc[~data.index.str.contains('Zernike'), :]
                data = data.set_index("Objects", drop=False)
                data.Mean = data.Mean.astype(float)
                data.Median = data.Median.astype(float)
                data.Std = data.Std.astype(float)

                if data.index.unique().__contains__('Cells'):
                    Cells.append(np.asarray(data.loc["Cells", "Mean"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cells", "Median"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cells", "Std"], dtype=np.float))
                elif data.index.unique().__contains__('Cell'):
                    Cells.append(np.asarray(data.loc["Cell", "Mean"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cell", "Median"], dtype=np.float))
                    Cells.append(np.asarray(data.loc["Cell", "Std"], dtype=np.float))

                Nuclei.append(np.asarray(data.loc["Nuclei", "Mean"], dtype=np.float))
                Nuclei.append(np.asarray(data.loc["Nuclei", "Median"], dtype=np.float))
                Nuclei.append(np.asarray(data.loc["Nuclei", "Std"], dtype=np.float))

                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Mean"], dtype=np.float))
                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Median"], dtype=np.float))
                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Std"], dtype=np.float))
            elif file.endswith('Image_MeasureTexture.csv'):
                #print(path+'/' + dir+'/'+file)
                sheet = pd.read_csv(path + '/' + dir + '/' + file, dtype=str)
                data = sheet.loc[:, ['Object', 'Value']]
                data = data.set_index("Object", drop=False)
                data.Value = data.Value.astype(float)

                if data.index.unique().__contains__('Cells'):
                    Cells.append(np.asarray(data.loc["Cells", "Value"], dtype=np.float))
                elif data.index.unique().__contains__('Cell'):
                    Cells.append(np.asarray(data.loc["Cell", "Value"], dtype=np.float))
                Nuclei.append(np.asarray(data.loc["Nuclei", "Value"], dtype=np.float))
                Cytoplasm.append(np.asarray(data.loc["Cytoplasm", "Value"], dtype=np.float))

        cells = np.array([elem for singleList in Cells for elem in singleList])
        nuclei = np.array([elem for singleList in Nuclei for elem in singleList])
        cytoplasm = np.array([elem for singleList in Cytoplasm for elem in singleList])

        all_cells.append(cells)
        all_nuclei.append(nuclei)
        all_cytoplasm.append(cytoplasm)
    all_cells = np.asarray(all_cells)
    all_nuclei = np.asarray(all_nuclei)
    all_cytoplasm = np.asarray(all_cytoplasm)
    return all_cells, all_nuclei, all_cytoplasm

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise IOError("Not enough inputs")
    path1 = sys.argv[1]
    path2 = sys.argv[2]
    print("Processing")
    cells = np.zeros((1, 1480), dtype=np.float) * np.nan
    nuclei = np.zeros((1, 1480), dtype=np.float) * np.nan
    cytoplasm = np.zeros((1, 1480), dtype=np.float) * np.nan

    cells_5 = np.zeros((40, 37), dtype=np.float) * np.nan
    nuclei_5 = np.zeros((40, 37), dtype=np.float) * np.nan
    cytoplasm_5 = np.zeros((40, 37), dtype=np.float) * np.nan


    cells_23_t, nuclei_23_t, cytoplasm_23_t = parse_image_measurement(path1)
    cells_5_t, nuclei_5_t, cytoplasm_5_t = parse_image_measurement(path2)



    for i in range(cells_5_t.shape[1]):
        cells[0, i] = np.corrcoef(cells_23_t[:, i], cells_5_t[:, i])[0, 1]
        nuclei[0, i] = np.corrcoef(nuclei_23_t[:, i], nuclei_5_t[:, i])[0, 1]
        cytoplasm[0, i] = np.corrcoef(cytoplasm_23_t[:, i], cytoplasm_5_t[:, i])[0, 1]
    cells = np.reshape(cells, (40, 37))
    nuclei = np.reshape(nuclei, (40, 37))
    cytoplasm = np.reshape(cytoplasm, (40, 37))

    fig, ax = plt.subplots()
    im, cbar = heatmap(cells, ax=ax, cmap="magma_r", grid_width=1)
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.title('Cell')
    plt.show()

    fig, ax = plt.subplots()
    im, cbar = heatmap(nuclei, ax=ax, cmap="magma_r", grid_width=1)
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.title('Nuclei')
    plt.show()

    fig, ax = plt.subplots()
    im, cbar = heatmap(cytoplasm, ax=ax, cmap="magma_r", grid_width=1)
    texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.title('Cytoplasm')
    plt.show()

    plt.figure(figsize=(15, 15))
    plt.rcParams.update({'font.size': 22})
    plt.grid()
    plt.hist([cells.ravel(), nuclei.ravel(), cytoplasm.ravel()], bins=25, alpha=1.0,
             label=['Cell', 'Nuclei', 'Cytoplasm'], edgecolor='black')

    #plt.hist(cells.ravel(), bins=25, alpha=1.0, edgecolor='black', label='Cell')
    #plt.hist(nuclei.ravel(), bins=25, alpha=0.7, edgecolor='black', label='Nuclei')
    #plt.hist(cytoplasm.ravel(), bins=25, alpha=0.5, edgecolor='black', label='Cytoplasm')
    plt.legend(loc='upper left')

    plt.show()

