function [ PyIsGood, cmdout ] = runPythonsetup(name)
% cf=pwd
% cd(['./' path])
%path='jw/INPUTS';
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%run abaqus must be in inp folder
% fid = fopen(['./runstatus.txt'], 'w+');
% trials = (dir(fullfile( '*.odb')));
% addpath(genpath(path));
%trialsForSTO =struct2cell(trialsForSTO );
% nTrials =length(trials);
% load('DISP_ROT.mat','disp_rot')
%for i= 1:nTrials
  %  if contains(trials(i).name, '-FT')~=1
        %name=[trials(i).name(1:end-4)];
        currentFolder = pwd;
        parentFolder = fileparts(currentFolder);
        addpath(genpath([parentFolder filesep 'RESULTS']));
        if exist([parentFolder filesep 'RESULTS_IMP2' filesep  filesep  name ], 'dir')==0
            mkdir([parentFolder filesep 'RESULTS_IMP2' filesep  filesep  name ])
        end
        
        if exist([parentFolder filesep 'RESULTS_IMP3' filesep  filesep  name ], 'dir')==0
            mkdir([parentFolder filesep 'RESULTS_IMP3' filesep  filesep  name ])
        end


%         end
%         try
         %   if contains(trials(i).name, 'initial')~=1
                
                pythonfile='GET_MOTIONS.py';
                odbfile=[name '.odb'];
                
                cmd2=['abaqus python ' pythonfile ' ' odbfile];
                [ PyIsGood, cmdout ]=system(cmd2);
                
                if PyIsGood==0
                   % fprintf(fid,['Trial: ' name ' complete\n']);
                    disp(['Trial: ' name ' complete'])
                else
                  %  fprintf(fid,['Trial: ' name ' failed\n']);
                    disp(cmdout);
                end
                
                file='MOTIONS.results';
                fid3 = fopen(file,'r');
                
                Z = textscan(fid3, '%s', 'delimiter', '\n'); % import full file
                
                % Z = textscan(fid, '\t%f\t%f\t%f\n','HeaderLines', 1); % import full file
                z = [Z{1}];
                fclose(fid3); % close file
                
                
                
                f1=fullfile([parentFolder filesep 'RESULTS_IMP1' filesep  name filesep],'target_AP.inp');
                f2=fullfile([parentFolder filesep 'RESULTS_IMP1' filesep  name filesep],'target_IE.inp');
                
                
                [fid1,errmsg1]  = fopen(f1,'w+');
                
                [fid2,errmsg2] = fopen(f2,'w+');
                formatSpec = '%4.3f\n';
               % first 0.25sec for stabilising
                initial=str2num(z{2});
                
                 for k=1:10
%                     fprintf(fid1,formatSpec,initial(2));
%                     fprintf(fid2,formatSpec,initial(3)*57.2958);
                       fprintf(fid1,formatSpec,0);
                    fprintf(fid2,formatSpec,0);
                end
               
                for j = 2:size(z,1)
                    C= str2num(z{j});
                    
%                     fprintf(fid1,formatSpec,C(2));
%                     fprintf(fid2,formatSpec,C(3))%*57.2958);
                     fprintf(fid1,formatSpec,C(2)-initial(2));
                    fprintf(fid2,formatSpec,(C(3)-initial(3)));%*57.2958);
                end
%                 for k=1:10
%                     fprintf(fid1,formatSpec,0);
%                     fprintf(fid2,formatSpec,0);
%                 end
%                
%                 for j = 2:size(z,1)
%                     C= str2num(z{j});
%                     
%                     fprintf(fid1,formatSpec,C(2)-initial(2));
%                     fprintf(fid2,formatSpec,C(3)-initial(3));
%                 end
%                % D=str2num(z{size(z,1)})-str2num(z{2});
%                 disp_rot.IE.(name)=initial(3);
%                    disp_rot.AP.(name)=initial(2);
                fclose(fid1);
                fclose(fid2);
                
                f3=fullfile([parentFolder filesep 'RESULTS_IMP2' filesep  name filesep],'target_AP.inp');
                f4=fullfile([parentFolder filesep 'RESULTS_IMP2' filesep  name filesep],'target_IE.inp');
                
                f5=fullfile([parentFolder filesep 'RESULTS_IMP3' filesep  name filesep],'target_AP.inp');
                f6=fullfile([parentFolder filesep 'RESULTS_IMP3' filesep  name filesep],'target_IE.inp');
                f7=fullfile([parentFolder filesep 'RESULTS_IMP4' filesep  name filesep],'target_AP.inp');
                f8=fullfile([parentFolder filesep 'RESULTS_IMP4' filesep  name filesep],'target_IE.inp');
                f9=fullfile([parentFolder filesep 'RESULTS_IMP5' filesep  name filesep],'target_AP.inp');
                f10=fullfile([parentFolder filesep 'RESULTS_IMP5' filesep  name filesep],'target_IE.inp');

                copyfile(f1,f3);
                copyfile(f2,f4);
                
                copyfile(f1,f5);
                copyfile(f2,f6);
                
                copyfile(f1,f7);
                copyfile(f2,f8);
                
                copyfile(f1,f9);
                copyfile(f2,f10);
                

                
        
                
                %%     DUMMY DIAG files
%                 f3=fullfile([parentFolder filesep 'RESULTS' filesep  name filesep],'diag-ie.txt');
%                 f4=fullfile([parentFolder filesep 'RESULTS' filesep  name filesep],'diag-ap.txt');
%                 
%                 if exist(f3,'file')==0
%                 [fid3,errmsg1]  = fopen(f3,'w+');
%                 
% %                 [fid4,errmsg2] = fopen(f4,'w+');
% %                 end
%                 
%                 fclose(fid3);
%                 fclose(fid4);
                
%             end
            
            
%         catch 
            %disp('error'); disp(errmsg1); disp(errmsg2);
            
%         end
        
        
%     end
   
% end
% save('DISP_ROT.mat','disp_rot')
% fclose(fid);
end


