function [ PyIsGood, cmdout ] = runPythonsetupinitial
% cf=pwd
% cd(['./' path])
%path='jw/INPUTS';
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%run abaqus must be in inp folder
cdir=pwd;
implants=1:6;
t={'ngait_og1'};% 'lungef1' 'stairdownf2' 'stairupf1' 'stepupf2'};
for i=1:length(implants);
    imp=implants(i);
    cd([cdir filesep 'IMP' num2str(imp)]);
    
   load( ['DISP_ROT_IMP' num2str(imp) '.mat'],'disp_rot')
fid = fopen(['./runstatus.txt'], 'w+');
trials = (dir(fullfile( '*.odb')));


% addpath(genpath(path));
%trialsForSTO =struct2cell(trialsForSTO );
nTrials =length(trials);
disp_rot=struct;
for i= 1:nTrials
  
        currentFolder = pwd;
        parentFolder = fileparts(currentFolder);
      
     
            if contains(trials(i).name, 'initial')==1
                  if contains(trials(i).name, '-FT')==1
        name=[trials(i).name(1:end-7)];
         odbfile=[name '-FT.odb'];
    else
        name=[trials(i).name(1:end-4)];
         odbfile=[name '.odb'];
    end
          %name=regexprep(name, '6', '7')   ;   
                pythonfile='GET_MOTIONSinitial.py';
               
                
                cmd2=['abaqus python ' pythonfile ' ' odbfile];
                [ PyIsGood, cmdout ]=system(cmd2);
                
                if PyIsGood==0
                    fprintf(fid,['Trial: ' name ' complete\n']);
                    disp(['Trial: ' name ' complete'])
                else
                    fprintf(fid,['Trial: ' name ' failed\n']);
                    disp(cmdout);
                end
                
                file='MOTIONS.results';
                fid3 = fopen(file,'r');
                
                Z = textscan(fid3, '%s', 'delimiter', '\n'); % import full file
                
                % Z = textscan(fid, '\t%f\t%f\t%f\n','HeaderLines', 1); % import full file
                z = [Z{1}];
                fclose(fid3); % close file
                
                
              
                initial=str2num(z{end});
                
              
                disp_rot.IE.(name)=initial(3);
                   disp_rot.AP.(name)=initial(2);
    

        end
        
        
  
   
end
save(['DISP_ROT_IMP' num2str(imp) '.mat'],'disp_rot','-append');
cd(cdir)
end
 fclose(fid);
end



