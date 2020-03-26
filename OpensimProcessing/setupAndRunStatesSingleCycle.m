function [  ] = setupAndRunStatesSingleCycle(varargin)
%Find states of the model over given time period. 
import org.opensim.modeling.*;
% move to directory where this subject's files are kept
patient=varargin{1};
path=varargin{2};
trial=[varargin{3} '_F'];
if length(varargin)==5
    initial_time=varargin{4};
    final_time=varargin{5};
    trial=[varargin{3}]; 
end

%model folder
model_folder= [path patient '/Scaling/'];

% specify where setup files will be printed.
state_setup_path = [path  patient '/StatesSetup/'];

% Model
scale_path = [path patient '/Scaling/'];
if strncmp(trial,'stair',5)==1 || strncmp(trial,'squat',5)==1 ...
        || strncmp(trial,'step',4)==1 || strncmp(trial,'chair',5)==1 || strncmp(trial,'lunge',5)==1
    
    model_file = [scale_path patient '_Model_fluro.osim'];
elseif  strncmp(trial,'dum',3)==1
    model_file = [scale_path patient '_Model_dummy.osim'];
else
    model_file = [scale_path patient '_Model.osim'];
end
model = Model(model_file);

%ik
ik_path =[path patient '/IKResults/'];
kin_file = [ik_path trial '_ik.mot']; 

% specify where results will be printed.
results_folder = [path  patient '/StatesResults/'];

name = [patient trial];
%%

    
    if length(varargin)==3 %get times
        motCoordsData=dlmread(kin_file,'\t',11,0);
        
        % for this example, column is time
        initial_time = motCoordsData(1,1);
        final_time=motCoordsData(end,1);
    end
    % get .mot data to determine time range

%Setup file
xmlDoc = xmlread([state_setup_path 'SetUp_SR_generic.xml']);
%change variable names
xmlDoc.getElementsByTagName('AnalyzeTool').item(0).getAttributes.item(0).setValue(trial);
xmlDoc.getElementsByTagName('model_file').item(0).getFirstChild.setNodeValue(model_file);
xmlDoc.getElementsByTagName('results_directory').item(0).getFirstChild.setNodeValue(results_folder);
xmlDoc.getElementsByTagName('coordinates_file').item(0).getFirstChild.setNodeValue(kin_file);
xmlDoc.getElementsByTagName('initial_time').item(0).getFirstChild.setNodeValue([num2str(initial_time)]);
xmlDoc.getElementsByTagName('final_time').item(0).getFirstChild.setNodeValue([num2str(final_time)]);
xmlDoc.getElementsByTagName('start_time').item(0).getFirstChild.setNodeValue([num2str(initial_time)]);
outfile = ['Setup_SR_' name '.xml'];
xmlwrite([state_setup_path '\' outfile], xmlDoc);


% Run SR
    

% newsrTool = AnalyzeTool([state_setup_path outfile]);
% newsrTool.run();


    
    %Execute
    Command = ['"c:\OpenSim 3.3\bin\analyze.exe" -S ' state_setup_path outfile];
    fprintf(['Performing states analysis on ' trial '\n']);
    [x,y]=system(Command); %stops printing
%     newsrTool = AnalyzeTool([setupfiles_folder outfile]);
% newsrTool.run();
    % rename the out.log so that it doesn't get overwritten
    copyfile('out.log',[results_folder '\' trial '_out.log'],'f');

end

