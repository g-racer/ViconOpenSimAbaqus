
%%
function setupAndRunIKSingleCycle(varargin)
%patient,path, trial, initial_time, final_time
import org.opensim.modeling.*;
% move to directory where this subject's files are kept
patient=varargin{1};
path=varargin{2};
trial=varargin{3};
if length(varargin)==5
    %trial=[varargin{3} '_SC']; %single cycle analysis
    initial_time=varargin{4};
    final_time=varargin{5};
end

%model folder
model_folder= [path patient '/Scaling/'];

% Go to the folder in the subject's folder where .trc files are
trc_data_folder = [path  patient '/MarkerData/'];
%uigetdir(subjectDir, 'Select the folder that contains the marker data files in .trc format.')

% specify where setup files will be printed.
setupfiles_folder = [path  patient '/IKSetup/'];

% specify where results will be printed.
results_folder = [path  patient '/IKResults/'];



%% Get and operate on the files
% choose a generic setup file to work from

    GenericSetupForIK = [setupfiles_folder 'Setup_IK_generic.xml'];
   xmlDoc = xmlread(GenericSetupForIK);

    %get the name of the file for this trial
    markerfile = [trc_data_folder trial '.trc'];
    
   
    % get trc data to determine time range
    % this example uses trc data that where column 0 is
    % frame number and column 1 is time.
    % verify that the row and column values (6 and 1 in this case) in
    % dlmread are appropriate for your input files.
    trcData=dlmread(markerfile,'\t',6,1);
       if length(varargin)==3
           % get initial and inal time from first and last entry in time column
    initial_time = trcData(1,1);
    final_time=trcData(end,1);
    trial=[trial '_F']; %full cycle analysis
       end
       
       if strncmp(trial,'stair',5)==1 || strncmp(trial,'squat',5)==1 ...
               || strncmp(trial,'step',4)==1 || strncmp(trial,'chair',5)==1 || strncmp(trial,'lunge',5)==1
           
           model_file = [model_folder patient '_Model_fluro.osim'];
       elseif  strncmp(trial,'dum',3)==1
           model_file = [model_folder patient '_Model_dummy.osim'];
       else
           model_file = [model_folder patient '_Model.osim'];
       end

    xmlDoc.getElementsByTagName('IKTool').item(0).getAttributes.item(0).setValue(trial);
    xmlDoc.getElementsByTagName('model_file').item(0).getFirstChild.setNodeValue(model_file);
    xmlDoc.getElementsByTagName('results_directory').item(0).getFirstChild.setNodeValue([results_folder]);
    xmlDoc.getElementsByTagName('marker_file').item(0).getFirstChild.setNodeValue(markerfile);
    xmlDoc.getElementsByTagName('output_motion_file').item(0).getFirstChild.setNodeValue([results_folder trial '_ik.mot']);
    xmlDoc.getElementsByTagName('time_range').item(0).getFirstChild.setNodeValue([num2str(initial_time) ' ' num2str(final_time)]);
    
    
    outfile = ([setupfiles_folder 'Setup_IK_' trial '.xml']);
    xmlwrite([outfile], xmlDoc);
    
%    newsrTool = IKTool(outfile);
%     newsrTool.run()
    Command = ['"c:\OpenSim 3.3\bin\ik.exe" -S ' outfile];
    fprintf(['Performing IK on cycle # ' num2str(trial) '\n']);
    [x, y]=system(Command)
    
    % rename the out.log so that it doesn't get overwritten
    copyfile('out.log',[results_folder trial '_out.log'],'f')
    
    
    
end
%sendmail(mail,subjectNumber, ['Hello! IK for ' subjectNumber '  is complete']);


