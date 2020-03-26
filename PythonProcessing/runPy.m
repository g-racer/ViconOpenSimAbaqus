function [pyIsGood, cmdout] = runPy(pyPath,odbPath,varargin)
%RUNPY Short description
%    [pyIsGood, cmdout] = runPy(pyPath,odbPath) will run the Python script
%    specified by pyPath on the Abaqus output database file specified by
%    odbPath. An output directory named PythonOutput will be produced in
%    the current working directory, and the result file will have the same
%    name as the Python script file, but a .result extension. Since the
%    python script is run through the system command line, a boolean value
%    pyIsGood will return a 1 if the script was successful, and the system
%    output from the command will be stored in cmdout.
%
%    Example:
%        [pyIsGood, cmdout] = runPy(pyPath,odbPath,outDirPath) will do the
%        same as above, but will place the .result files in the directory
%        specified by outDirPath. If the directory does not exist, it will
%        be created.
%
%    Note: 
%        < Removed - See source comments.>
%        If the python script specified by pyPath is GET_COORDS.py, runPy
%        will convert the coordinate system to a Grood-Suntay joint and
%        output kinematic results in addition to the normal script output.
%
%    See also SYSTEM, BUILDGSSTRUCT, GROODSUNTAYJOINT, SAVEKINEMATICRESULTS 

%
%---------------------------------------------------------------------------
p = inputParser;

% Define parser defaults.
defaultPyOutDir = fullfile('.','PythonOutput');

% Add variable definitions to input parser object.
addRequired(p,'pyPath',@ischar);
addRequired(p,'odbPath',@ischar);
addOptional(p,'pyOutDir',defaultPyOutDir,@ischar);

% Parse the inputs and collect the variables.
parse(p,pyPath,odbPath,varargin{:});
pyOutDir = p.Results.pyOutDir;
pyPath = p.Results.pyPath;
odbPath = p.Results.odbPath;


%% Make sure the output directory exists
if ~exist(pyOutDir,'dir')
    mkdir(pyOutDir);
    addpath(genpath(pyOutDir));
end


%% Break out the filenames from their paths & build other names
[~,odbFileName,odbExt] = fileparts(odbPath);
[~,pyFileName,~] = fileparts(pyPath);
pyResultName = [pyFileName, '.results'];
pyResultPath = fullfile(pyOutDir,pyResultName);


%% Set up the system command to run the script
command = [...
           'abaqus python "',...
           pyPath,...
           '" ',...
           [odbFileName,odbExt],...
           ' ',...
           pyResultPath...
          ];
      
      
%% Run the python script through terminal.
[pyHasFailed, cmdout] = system(command);

% Return an affirmitive bool name
pyIsGood = ~pyHasFailed; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This piece of code was specific to Kalin's project. Ask about it if you 
%% are just starting to set up Grood-Suntay joints in Matlab.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
% Handle the GET_COORDS file, making kinematic results.
if contains(pyResultPath,'GET_COORDS') && pyIsGood
    % Build the Grood-Suntay structure.
    GS = buildGsStruct(...
        pyResultPath,...
        'right',...
        {'Fem','Tib','Pat'}...
    );
    
    % Set up the Grood-Suntay joints.
    Tibia = groodSuntayJoint(GS,'Fem','Tib');
    Patella = groodSuntayJoint(GS,'Fem','Pat');
    saveKinematicResults(Tibia,pyOutDir);
    saveKinematicResults(Patella,pyOutDir);
end
%}

end