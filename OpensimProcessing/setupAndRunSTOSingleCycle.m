%run STO on a single trial

function setupAndRunSTOSingleCycle(varargin)
 import org.opensim.modeling.* % % Import OpenSim Libraries
% move to directory where this subject's files are kept

patient=varargin{1};
path=varargin{2};
trial=[varargin{3} '_F'];
if length(varargin)==5
    initial_time=varargin{4};
    final_time=varargin{5};
    trial=[varargin{3}]; 
end

% Folder in the subject's folder where .trc IK result files are
ik_results_folder = [path patient '/IKResults/'];
%uigetdir(subjectDir, 'Select the folder that contains the marker data files in .trc format.')

% folder containing force plate data
Forcedata_folder=[path patient '/ForceData/'];

% specify where setup files will be printed.
setupfiles_folder = [path patient '/STOSetup/'];
%uigetdir(subjectDir, 'Select the folder where the IK Setup Files will be printed.');

% specify where results will be printed.
results_folder = [path patient '/STOresults/'];
%uigetdir(subjectDir, 'Select the folder where the IK Results will be printed.')

% specify where EL setup files will be printed.
ELsetupfiles_folder=[path patient '/ELSetup/'];

%model folder
model_folder= [path patient '/Scaling/'];
 if strncmp(trial,'stair',5)==1 || strncmp(trial,'squat',5)==1 ...
            || strncmp(trial,'step',4)==1 || strncmp(trial,'chair',5)==1 || strncmp(trial,'lunge',5)==1
        
model_file = [model_folder patient '_Model_fluro.osim'];
else
    model_file = [model_folder patient '_Model.osim'];
       end




%% Get and operate on the files
% choose files in results folder with .csv ending

genericSetupForAn = [setupfiles_folder 'Setup_Analyze_generic.xml'];
xmlDoc = xmlread(genericSetupForAn);
xmlDoc.getElementsByTagName('model_file').item(0).getFirstChild.setNodeValue(model_file);

% get the file names that match the ik_reults convention
% this is where consistent naming conventions pay off
   
    %single cycle analysis
    
    ikfile=[ik_results_folder trial '_ik.mot'];
    
    if length(varargin)==3 %get times
        motCoordsData=dlmread([ikfile],'\t',11,0);
      
        initial_time = motCoordsData(1,1);
        final_time=motCoordsData(end,1);
    end
    % get .mot data to determine time range

     xmlDoc.getElementsByTagName('AnalyzeTool').item(0).getAttributes.item(0).setValue(trial);
    xmlDoc.getElementsByTagName('model_file').item(0).getFirstChild.setNodeValue([model_file]);
    xmlDoc.getElementsByTagName('results_directory').item(0).getFirstChild.setNodeValue(results_folder);
    xmlDoc.getElementsByTagName('coordinates_file').item(0).getFirstChild.setNodeValue(ikfile);
    xmlDoc.getElementsByTagName('initial_time').item(0).getFirstChild.setNodeValue([num2str(initial_time)]);
    xmlDoc.getElementsByTagName('final_time').item(0).getFirstChild.setNodeValue([num2str(final_time)]);
    xmlDoc.getElementsByTagName('start_time').item(0).getFirstChild.setNodeValue([num2str(initial_time)]);
    xmlDoc.getElementsByTagName('end_time').item(0).getFirstChild.setNodeValue([num2str(final_time)]);
    
     ELoutfile = ['EL_Setup_' trial '.xml'];
        %writes EL setup file
           
    
    %if strncmp(trial,'L',1) || strncmp(trial,'R',1)==1 %if gait trial with force plate data
        if strncmp(trial,'RWalk',5)
            ELSetupForSTO = [ELsetupfiles_folder patient '_RWalk_ELSetup.xml'];
        elseif strncmp(trial,'RCut',4) || strncmp(trial,'RPivot',5)
            ELSetupForSTO = [ELsetupfiles_folder 'R_ELSetup.xml'];
        elseif strncmp(trial,'LCut',4) || strncmp(trial,'LPivot',5)
            ELSetupForSTO = [ELsetupfiles_folder 'L_ELSetup.xml'];
            
            
%             elseif strncmp(trial,'stairup',7)==1 
%                  ELSetupForSTO = [ELsetupfiles_folder 'R_12_ELSetup.xml'];
                 
                elseif strncmp(trial,'stairdo',7)==1 || strncmp(trial,'stairup',7)==1  ||  ...
                        strncmp(trial,'lunge',5)==1 || strncmp(trial,'chair',5)==1 ...
                   || strncmp(trial,'step',4)==1
                                 ELSetupForSTO = [ELsetupfiles_folder 'R_1_ELSetup.xml'];  
                 
                  elseif   strncmp(trial,'chair',5)==1 ...  
                        || strncmp(trial,'twist',5)==1 
                 ELSetupForSTO = [ELsetupfiles_folder 'R_3_ELSetup.xml'];
                 
%                   elseif  strncmp(trial,'step',4)==1 
%                  ELSetupForSTO = [];
                 
        elseif contains(trial,'gait')==1 || strncmp(trial,'bouncy',6)==1 ...
            || contains(trial,'crouch')==1 || strncmp(trial,'med',3)==1
        ELSetupForSTO = [ELsetupfiles_folder 'R_gait_ELSetup.xml']; %FP 2 and 3
        end
        %setup external loads file with correct force data from specific trial
        xmlDoc_EL = xmlread(ELSetupForSTO);
        xmlDoc_EL.getElementsByTagName('datafile').item(0).getFirstChild.setNodeValue([Forcedata_folder trial '.mot']);
        
        %Writes output setup file for external loads
        xmlwrite([ELsetupfiles_folder '/' ELoutfile], xmlDoc_EL)
        
        xmlDoc.getElementsByTagName('external_loads_file').item(0).getFirstChild.setNodeValue([ELsetupfiles_folder  ELoutfile]);
    %end
        
    outfile = ['Setup_Analyze_' trial '.xml'];
    xmlwrite([setupfiles_folder outfile], xmlDoc);
    %Execute
    Command = ['"c:\OpenSim 3.3\bin\analyze.exe" -S ' setupfiles_folder outfile];
    fprintf(['Performing STO on ' trial '\n']);
    [x,y]=system(Command); %stops printing
%     newsrTool = AnalyzeTool([setupfiles_folder outfile]);
% newsrTool.run();
    % rename the out.log so that it doesn't get overwritten
    copyfile('out.log',[results_folder '/' trial '_out.log'],'f');
    
    
end   

%sendmail(mail,subjectNumber, ['Hello! Analysis for ' subjectNumber '  is complete']);


