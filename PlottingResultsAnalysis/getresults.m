function [ RESULTS ] = getresults
types={'lunge', 'gait', 'staird', 'stairup', 'stepup'};
nos=[6]%1 12 13 14 10 ];

original_dir=pwd;
load('results.mat','RESULTS')
for t=1:length(types)
    type=types{t};
    for o=1:length(nos)
        num=num2str(nos(o));
        
        path=['jw' filesep type num];
        
        cd(path)
        trials = (dir(fullfile('*.odb')));
        
        nTrials =length(trials);
        
        
        for i= 1:nTrials
            
            currentFolder = path; 
            if contains(trials(i).name, 'initial')==0
                if contains(trials(i).name, '-FT')==1
                    name=[trials(i).name(1:end-7)];
                    odbfile=[name '-FT.odb'];
                else
                    name=[trials(i).name(1:end-4)];
                    odbfile=[name '.odb'];
                end
                
                ind=strfind(name,'step');
                ind=ind(end);
                trial=name(1:ind-1);
                step=name(ind+4);
                ld=name(ind+6);
                try
                    implant=(name(ind+11));
                catch
                    implant='0';
                end
                            
                pythonfiles=dir(fullfile( '*.py'));
                Py={pythonfiles(:).name};
                Py(contains(Py(:, :), 'GET_FORCE_MUSC')) = [];
                Py(contains(Py(:, :), 'GET_KIN_KNEE')) = [];
                Py(contains(Py(:, :), 'GET_JOINT')) = [];
                Py{end+1}=['GET_FORCE_MUSC' ld '.py'];
                
                for k=1:length(Py)
                    pythonfile=Py{k};
                    
                    p=strfind(pythonfile,'.py');
                    resulttype=pythonfile(5:p-1);
                    if contains(resulttype,'PATIENT')
                        resulttype=    erase(resulttype,"_PATIENT");
                    end
                    cmd2=['abaqus python ' pythonfile ' ' odbfile];
                    [ PyIsGood, cmdout ]=system(cmd2);
                    
                    if PyIsGood==0
                        % fprintf(fid,['Trial: ' name ' complete\n']);
                        disp(['Trial: ' name ' ' resulttype ' script complete'])
                    else
                        %fprintf(fid,['Trial: ' name ' failed\n']);
                        disp(cmdout);
                    end
                    
                    
                    
                    file=[resulttype '.results'];
                    
                    try
                        Data=dlmread(file,'',1,0);
                        
                        %Data=D.data;
                        if isempty(Data)==1
                            disp([name ' is empty'])
                            continue
                        else
                            if contains(resulttype, 'FORCE_MUSC')==1
                                resulttype='FORCE_MUSC';
                            end
                            
                            RESULTS.(resulttype).(trial).(['step' step]).(ld).(['implant' implant])=Data;
                            % RESULTS.(resulttype).('headers')=D.colheaders;
                        end
                        
                    catch
                    end
                    
                    
                end
                disp([name ' completed'])
            end
            
            
        end
        
        cd(original_dir)
        
        
        
        
    end
    
end

save([original_dir filesep  'results.mat'],'RESULTS','-append');
end




