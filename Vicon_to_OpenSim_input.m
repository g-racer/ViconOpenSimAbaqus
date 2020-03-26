%Setup file to run the Vicon to OpenSim MATLAB function.
clear all
close all
path= 'C:\Users\Gracemcconnochie\Google_Drive\Thesis\ACL_Pilot_data\Processing\';
patient='VA' ;
csv_folder = [path patient '\CSV_files\'];
markerdata_folder = [path patient '\MarkerData\'];
forcedata_folder= [path patient '\ForceData\'];


trialsForAn = dir([csv_folder '*.csv']);
nTrials =length(trialsForAn);

direction = 1;
vicon_version = '2.3';
cd(path)
for trial= 1:nTrials
    %get the name of the file for this trial
    csv_trial = trialsForAn(1).name;
    
    % create name of trial from .csv file name
    trial = regexprep(csv_trial,'.csv','');
    csv_file = [csv_folder csv_trial] ;
    mot_file = [forcedata_folder trial '.mot'] ;
    trc_file = [markerdata_folder trial '.trc'] ;
    %Run the function
    Vicon_to_OpenSim(csv_file,trc_file,mot_file,vicon_version,direction)
end


