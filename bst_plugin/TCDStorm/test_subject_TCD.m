close all
clear all
clc

fname='RC071_T6_rest_stand_TCD_RMCA_LMCA'; %participant file name without .mat
path_to_LabChart_data='C:\Users\ben82\OneDrive - Universite de Montreal\Stage-H2024\Recardio_Analysis_April_2024\Recardio_TCD_April2024\RC_TCD_Data';
TCD_manifest_folder_path = 'C:\Users\ben82\OneDrive - Universite de Montreal\Stage-H2024\Recardio_Analysis_April_2024\Recardio_TCD_April2024\RC_TCD_manifest';
FingerBP_result_folder_path = 'C:\Users\ben82\OneDrive - Universite de Montreal\Stage-H2024\Recardio_Analysis_April_2024\Recardio_TCD_April2024\RC_FingerBP_results'; % Modify this path as necessary
MCA_result_folder_path ='C:\Users\ben82\OneDrive - Universite de Montreal\Stage-H2024\Recardio_Analysis_April_2024\Recardio_TCD_April2024\RC_TCD_results';
AutoReg_result_folder_path='C:\Users\ben82\OneDrive - Universite de Montreal\Stage-H2024\Recardio_Analysis_April_2024\Recardio_TCD_April2024\RC_AutoReg_results';

%% this analysis aussumes that TCD manifest is filled put but  findpeak values are defult it can be changed if its not good

GetDataInfo_From_assignedName(TCD_manifest_folder_path, fname)
fetch_TCD_data(path_to_LabChart_data, TCD_manifest_folder_path, fname);
verify_TCDdata_units(TCD_manifest_folder_path,fname);
ECG_RwaveLocations(TCD_manifest_folder_path,fname);
analyse_FingerBP(TCD_manifest_folder_path,fname, FingerBP_result_folder_path)
analyse_LMCA_RMCA(TCD_manifest_folder_path,fname, MCA_result_folder_path)
analyse_autoregulation(TCD_manifest_folder_path, fname, AutoReg_result_folder_path,0) % 10 is the 10-second movi average window
%calc_indicies_TCD_BP(signal, RwaveLocations, signalName)
%% TO DO: manual correction of ECG should have a function




 