import os
import numpy as np
import pandas as pd
import glob
import defaults as I
import sys

def MakeCorrelationCsv(working_dir = None, out_dir = None, mode = 1, measureStd = False):
    if working_dir is None:
        working_dir = I.working_dir
    if out_dir is None:
        out_dir = working_dir
    object_list = ["Nuclei.csv", "Cells.csv", "Cytoplasm.csv"]
    prefix_file_mod = 'Correlation_Correlation'
    dataFrame = []
    for pair in I.list_of_correlation_pairs:

        for obj in object_list[:3]:
            df_groupObject= (pd.read_csv(os.path.join(working_dir, I.prefix_csv_name + obj)))
            # print(I.prefix_csv_name + obj)
            # print(pair)
            col_name = prefix_file_mod + "_" + pair[0]+ "_" +  pair[1]
            try:
                # A = [i for i in df_groupObject.keys() if (pair[0] in i) and (pair[1] in i) and prefix_file_mod in i]
                col_data_temp = df_groupObject[col_name]
            except:
                ValueError("Print no value")
                return 0
            if mode == 1:
                temp_list = [pair[0], pair[1], obj.split('.')[0],
                             round(col_data_temp.mean(),2),
                             round(col_data_temp.median(),2),
                             round(col_data_temp.min(),2),
                             round(col_data_temp.max(),2),
                             round(col_data_temp.std(ddof = 0),2)]
                dataFrame.append(temp_list)
            else:
                dataFrame.append([pair[0], pair[1], obj.split('.')[0], "Mean Correlation",
                             round(col_data_temp.mean(),2)])
                dataFrame.append([pair[0], pair[1], obj.split('.')[0], "Median Correlation",
                                  round(col_data_temp.median(), 2)])
                dataFrame.append([pair[0], pair[1], obj.split('.')[0], "Min Correlation",
                                  round(col_data_temp.min(), 2)])
                dataFrame.append([pair[0], pair[1], obj.split('.')[0], "Max Correlation",
                                  round(col_data_temp.max(), 2)])

                # dataFrame.append([pair[0], pair[1], obj.split('.')[0], "Std",
                #                   round(col_data_temp.std(), 2)])
    if mode == 1:
        column_names = ['FirstImage', 'SecondImage','Object', 'Mean', 'Median', 'Min', 'Max','Std']
    else:
        column_names = ['FirstImage', 'SecondImage','Object', 'Measurements', 'Value']
    df = pd.DataFrame(dataFrame, columns=column_names)
    df.to_csv(os.path.join(out_dir,"Image_MeasureCorrelation.csv".format(mode)), index= False)
    return dataFrame


def MakeGranularityCsv(working_dir = None, out_dir = None):
    if working_dir is None:
        working_dir = I.working_dir
    if out_dir is None:
        out_dir = working_dir
    protein_list = ["Hoechst", "Mito","ERSyto", "PhGolgi", "ERSytoBleed"]
    number_of_granularity = 16
    columns_name = ["GS" + str(i) for i in range(1, number_of_granularity + 1)]
    columns_name = ["Protein"] + columns_name
    gran_prefix_name = "Granularity_"

    df_groupImage = pd.read_csv(os.path.join(working_dir, I.prefix_csv_name + "Image.csv"))
    frame = []

    for j in protein_list:
        gran_temp = []
        gran_temp.append(j)
        for i in range(1, number_of_granularity + 1):
            col_name = gran_prefix_name + str(i) + "_" + j
            col_data = df_groupImage[col_name]
            #         print(col_data.get_values)
            gran_temp.append(round(col_data[0], 2))
        frame.append(gran_temp)

    Measure_Granularity_files = pd.DataFrame(frame, columns=columns_name)
    Measure_Granularity_files.to_csv(os.path.join(out_dir,"Image_MeasureGranularity.csv".format(1)),
                                     index=False)
    return frame

def MakeObjectIntensityCsv(working_dir = None, out_dir = None, ddof= 0):
    if working_dir is None:
        working_dir = I.working_dir
    if out_dir is None:
        out_dir = working_dir

    object_list = ['Nuclei.csv', 'Cytoplasm.csv', 'Cells.csv', 'Image.csv']
    protein_list = ["Hoechst",
                    "ERSyto",
                    "ERSytoBleed",
                    "PhGolgi", "Mito"]
    frame = []
    for k in protein_list:
        for obj in object_list[:3]:
            groupObject = pd.read_csv(os.path.join(working_dir,
                                                   I.prefix_csv_name + obj))
            Obj = obj.split('.')[0]
            #         print(k, Obj)
            for j in I.ObjectIntensity_list:
                if ("_X" not in j) and ("_Y" not in j):
                    #                 print("\t",j)
                    col_data = groupObject["Intensity_" + j + "_" + k]
                    mean_col_data = round(col_data.mean(), 3)
                    median_col_data = round(col_data.median(), 3)
                    std_col_data = round(col_data.std(ddof = ddof), 3)

                    temp_series = [k, Obj, j,
                                   mean_col_data, median_col_data,
                                   std_col_data]
                    frame.append(temp_series)
                else:
                    col_data = groupObject["Location_" + j + "_" + k]
                    mean_col_data = round(col_data.mean(), 3)
                    median_col_data = round(col_data.median(), 3)
                    std_col_data = round(col_data.std(ddof =0), 3)
                    temp_series = [k, Obj, j,
                                   mean_col_data, median_col_data,
                                   std_col_data]
                    frame.append(temp_series)

    ObjectIntensity_file = pd.DataFrame(frame, columns=["Image", "Objects", "Feature",
                                                             "Mean", "Median",
                                                             "STD"], )
    ObjectIntensity_file.to_csv(os.path.join(out_dir,
                                                  "Image_MeasureObjectIntensity.csv"), index=False)
    return frame

def MakeObjectSizeShape(working_dir = None, out_dir = None, ddof = 0):
    if working_dir is None:
        working_dir = I.working_dir
    if out_dir is None:
        out_dir = working_dir

    frame = []
    for i in I.object_list[:3]:
        groupObject = pd.read_csv(os.path.join(working_dir,
                                               I.prefix_csv_name + i))
        Obj = i.split('.')[0]

        for j in I.ObjectSizeShape_list:
            for k in groupObject.keys():
                if ("AreaShape_" + j in k):
                    data_col = groupObject[k]
                    mean_data_col = round(data_col.mean(), 2)
                    median_data_col = round(data_col.median(), 2)
                    std_data_col = round(data_col.std(ddof = ddof), 2)

                    temp_series = [Obj, j, mean_data_col, median_data_col, std_data_col]
                    # print(temp_series)
                    frame.append(temp_series)

    Measure_SizeShape_files = pd.DataFrame(frame, columns=["Objects", "Feature",
                                                           "Mean", "Median",
                                                           "Std"], )
    Measure_SizeShape_files.to_csv(os.path.join(out_dir,
                                                "Image_MeasureObjectSizeShape"'.csv'), index=False)

def MakeObjectRadialDistributionCsv(working_dir = None, out_dir = None, mode = 2):
    if working_dir is None:
        working_dir = I.working_dir
    if out_dir is None:
        out_dir = working_dir
    protein_list = ["ERSyto",
                    "ERSytoBleed",
                    "PhGolgi", "Mito"]
    features = ['AtD', 'MeanFrac', 'RadialCV']
    frame = []
    if mode == 1:
        for i in I.object_list[:3]:
            groupObject = pd.read_csv(os.path.join(working_dir,
                                                   I.prefix_csv_name + i))
            for j in protein_list:
                count = 0
                for f in features:
                    if 'AtD' in f:
                        feat = 'Fraction'
                    elif 'MeanFrac' in f:
                        feat = 'Intensity'
                    elif 'RadialCV' in f:
                        feat = 'COV'
                    temp_series = [i.split('.')[0], j, feat]
                    for k in groupObject.keys():
                        if ("RadialDistribution" in k) and (j + "_" in k) and (f in k):
                            count += 1
                            data_col = groupObject[k]
                            temp_series.append(round(data_col.mean(), 4))
                        else:
                            continue
                    # print(temp_series)
                    frame.append(temp_series)
        Measure_ObjectRadialDistribution_files = pd.DataFrame(frame, columns=["Objects", "Protein", "Feature",
                                                                              "1to4", "2to4",
                                                                              "3to4", "4to4"])
        Measure_ObjectRadialDistribution_files.to_csv(os.path.join(out_dir,
                                                                   "Image_MeasureObjectRadialDistribution" + '.csv'),
                                                      index=False)
    else:
        for j in protein_list:
            for i in I.object_list[:3]:
                groupObject = pd.read_csv(os.path.join(working_dir,
                                                       I.prefix_csv_name + i))
                for bin in range(1,5):
                    AtD_key = "RadialDistribution_FracAtD_{}_{}of4".format(j,bin)
                    MeanFrac_key = "RadialDistribution_MeanFrac_{}_{}of4".format(j,bin)
                    Radial_key = "RadialDistribution_RadialCV_{}_{}of4".format(j,bin)
                    data_atd = round(groupObject[AtD_key].mean(),4)
                    data_mf = round(groupObject[MeanFrac_key].mean(),4)
                    temp = groupObject[Radial_key]
                    N = len(temp[temp>0])
                    data_rad = round(temp.sum()/N,4)
                    temp_series = [j,i.split('.')[0],  bin, data_atd, data_mf, data_rad]
                    frame.append(temp_series)


        Measure_ObjectRadialDistribution_files = pd.DataFrame(frame, columns=["Image", "Object", "bin",
                                                                              "Fraction", "Intensity",
                                                                              "COV"])
        Measure_ObjectRadialDistribution_files.to_csv(os.path.join(out_dir,
                                                                   "Image_MeasureObjectRadialDistribution" + '.csv'),
                                                      index=False)
    return frame

def MakeTextureCsv(working_dir = None, out_dir = None):
    if working_dir is None:
        working_dir = I.working_dir
    if out_dir is None:
        out_dir = working_dir
    sizes = [3, 5, 10]
    groupImage = pd.read_csv(os.path.join(working_dir,
                                          I.prefix_csv_name + 'Image.csv'))
    object_list = ['Cells.csv', 'Cytoplasm.csv', 'Nuclei.csv', 'Image.csv']
    protein_list = ["Hoechst",
                    "ERSyto",
                    "ERSytoBleed",
                    "PhGolgi", "Mito"]
    frame = []
    for pro in protein_list:
        for sz in sizes:
            Obj = None
            for feature in I.Texture_main_listA:
                name_key = "Texture_" + feature + "_" + pro + "_" + str(sz)
                # print(name_key)
                for j in groupImage.keys():
                    if name_key in j:
                        data_value = round(groupImage[j][0], 2)
                        temp = [pro, Obj, feature, sz, data_value]
                        frame.append(temp)

            for obj in object_list[:3]:
                groupObject = pd.read_csv(os.path.join(working_dir,
                                                       I.prefix_csv_name + obj))
                Obj = obj.split('.')[0]

                for feature in I.Texture_main_listB:
                    name_key = "Texture_" + feature + "_" + pro + "_" + str(sz)
                    for j in groupObject.keys():
                        if name_key in j:
                            data_value = groupObject[j]
                            temp = [pro, Obj, 'min ' +feature, sz, round(data_value.min(), 2)]
                            frame.append(temp)
                            temp = [pro, Obj, 'max ' +feature, sz, round(data_value.max(), 2)]
                            frame.append(temp)
                            temp = [pro, Obj, 'mean ' +feature, sz, round(data_value.mean(), 2)]
                            frame.append(temp)
                            temp = [pro, Obj, 'median ' +feature, sz, round(data_value.median(), 2)]
                            frame.append(temp)
                            temp = [pro, Obj, 'std dev '+ feature, sz, round(data_value.std(), 2)]
                            frame.append(temp)

    Measure_Texture_files = pd.DataFrame(frame, columns=["Image", "Object",
                                                         "Measurement", "Scale",
                                                         "Value"], )
    Measure_Texture_files.to_csv(os.path.join(out_dir,
                                              "Image_MeasureTexture.csv"), index=False)

    return frame
def get_list_of_image_a_ids(N = 100):
    ids_list = []
    gt_id_list = []
    a = 1
    count = 1
    for i in range(1, N + 1):
        gt_id = "Image{:03d}".format(i)
        ids_list.append([gt_id,"_a{:02d}_s{}".format(a, count)])
        count += 1
        if count == 10:
            count = 1
            a += 1
    return ids_list

def process_group_of_measurements(working_dir, output_dir):
    print("Measuring Correlation 1/6")
    MakeCorrelationCsv(working_dir, output_dir, mode=2)
    print("Measuring Granularity 2/6")
    MakeGranularityCsv(working_dir, output_dir)
    print("Measuring ObjectIntensity 3/6")
    MakeObjectIntensityCsv(working_dir, output_dir)
    print("Measuring RadialDistribution 4/6")
    MakeObjectRadialDistributionCsv(working_dir, output_dir)
    print("Measuring Object SizeShape 5/6")
    MakeObjectSizeShape(working_dir, output_dir)
    print("Measuring Object Texture 6/6")
    MakeTextureCsv(working_dir, output_dir)
    print(working_dir, "is Done")

def stop_when_not_found (path_, folder):
    if not os.path.exists(os.path.join(path_, folder)):
        raise IOError ("Not found")

def get_the_list_of_files_from_idList(id_list, input_folder):
    list_of_files = []
    for i in os.listdir(input_folder):
        for j in range(len(id_list)):
            if id_list[j][1] in i:
                list_of_files.append(os.path.join(input_folder, i))
    return list_of_files

def main():
    # The folder has all of the processed folders from Cell Profilers!
    working_dir = "D:\\100_protein_report\ouput_38034_original_data_CPv3_auto\ouput_auto_CPv3_38034_a01_s3"
    output_folder ="D:\\100_protein_report\\38034_original_data_process_CPv3"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get the list of "a01_s1, a01_s2, ... a01_s9,... a03_s9.... Default 100 files
    ids_list = get_list_of_image_a_ids()

    # Get all the folder that have matched IDs
    list_of_inputs = get_the_list_of_files_from_idList(ids_list, working_dir)
    if len(list_of_inputs) > 0:
        for input_dir in list_of_inputs:
            output_fn = input_dir.split("\\")[-1] + "_Auto"
            output_dir = os.path.join(output_folder, output_fn)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            process_group_of_measurements(input_dir, output_dir)
    else:
        input_dir = working_dir
        output_fn = input_dir.split("\\")[-1] + "_Auto"
        output_dir = os.path.join(output_folder, output_fn)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        process_group_of_measurements(input_dir, output_dir)


#main()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise IOError("Not enough inputs")

    WorkingDir = sys.argv[1]
    OutputDir = sys.argv[2]
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    ids_list = get_list_of_image_a_ids()
    list_of_inputs = get_the_list_of_files_from_idList(ids_list, WorkingDir)

    for input_dir in list_of_inputs:
        output_fn = input_dir.split("\\")[-1] + "_Auto"
        output_dir = os.path.join(OutputDir, output_fn)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        process_group_of_measurements(input_dir, output_dir)


