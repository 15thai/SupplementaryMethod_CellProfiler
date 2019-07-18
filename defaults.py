import os
import numpy as np
import pandas as pd
import glob

#working_dir = "D:\\100_protein_report\output_anh_a01_predict_test\\thanh_folder2_38034_a01_s1"
#prefix_folder_name = "thanh_folder2_38034"
#prefix_file_name = "MyExp_next_group"
#number_of_files = 100
prefix_csv_name = "MyExpt_next_group"
object_list = ['Cells.csv','Nuclei.csv','Cytoplasm.csv', 'Image.csv']
measure_list = ["IdentifyPrimary",
                "IdentifySecondary",
                "ObjectIntensity",
                "ObjectIntensityDistribution",
                "Texture",
                "SizeShape",
                "Colocalization",
                "Granularity"]

protein_list = ["Hoechst",
               "Mito",
               "ERSyto",
               "ERSytoBleed",
               "PhGolgi"]

list_of_correlation_pairs =[("Hoechst", "Mito"),
               ("Hoechst", "ERSyto"),
                ("Hoechst", "ERSytoBleed"),
                ("Hoechst","PhGolgi"),
               ("Mito","ERSyto"),
               ("Mito","ERSytoBleed"),
               ("Mito","PhGolgi"),
               ("ERSyto","ERSytoBleed"),
               ("ERSyto","PhGolgi"),
               ("ERSytoBleed","PhGolgi")]

correlation_measurement_columns= ["Objects", "First Image","Second Image",
                                  "Mean Correlation", "Median Correlation",
                                  "Min Correlation", "Max Correlation"]
ObjectIntensity_list = ["IntegratedIntensity",
                        "MeanIntensity",
                        "StdIntensity",
                        "MinIntensity",
                        "MaxIntensity",
                        "IntegratedIntensityEdge",
                        "MeanIntensityEdge",
                        "StdIntensityEdge",
                        "MinIntensityEdge",
                        "MaxIntensityEdge",
                        "MassDisplacement",
                        "LowerQuartileIntensity",
                        "MedianIntensity",
                        "MADIntensity",
                        "UpperQuartileIntensity",
                        "CenterMassIntensity_X",
                        "CenterMassIntensity_Y",
                        "MaxIntensity_X",
                        "MaxIntensity_Y"]

ObjectSizeShape_list = ["Eccentricity","MajorAxisLength","MinorAxisLength","Orientation","Compactness","Area","Center_X",
                   "Center_Y","Extent","Perimeter","Solidity","FormFactor","EulerNumber","MaximumRadius","MeanRadius",
                   "MedianRadius","MinFeretDiameter","MaxFeretDiameter","Zernike_0_0","Zernike_1_1","Zernike_2_0",
                   "Zernike_2_2","Zernike_3_1","Zernike_3_3","Zernike_4_0","Zernike_4_2","Zernike_4_4","Zernike_5_1",
                   "Zernike_5_3","Zernike_5_5","Zernike_6_0","Zernike_6_2","Zernike_6_4","Zernike_6_6","Zernike_7_1",
                   "Zernike_7_3","Zernike_7_5","Zernike_7_7","Zernike_8_0","Zernike_8_2","Zernike_8_4","Zernike_8_6",
                   "Zernike_8_8","Zernike_9_1","Zernike_9_3","Zernike_9_5","Zernike_9_7","Zernike_9_9"]

Texture_size_list = ["Gabor", "AngularSecondMoment", "Contrast", "Correlation", "Variance", "InverseDifferenceMoment",
                     "SumAverage","SumVariance","SumEntropy","Entropy","DifferenceVariance","DifferenceEntropy","InfoMeas1",
                     "InfoMeas2"]

Texture_main_listA = ["Gabor", "AngularSecondMoment", "Contrast", "Correlation", "Variance",
                      "InverseDifferenceMoment",
                      "SumAverage", "SumVariance", "SumEntropy", "Entropy", "DifferenceVariance",
                      "DifferenceEntropy", "InfoMeas1",
                      "InfoMeas2"]
Texture_main_listB = ["AngularSecondMoment", "Contrast", "Correlation", "Variance", "InverseDifferenceMoment",
                      "SumAverage", "SumVariance", "SumEntropy", "Entropy", "DifferenceVariance",
                      "DifferenceEntropy", "InfoMeas1",
                      "InfoMeas2", "Gabor", ]
