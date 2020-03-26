%path='./jw/RESULTS';
currentFolder = pwd;
path=[currentFolder  filesep 'jw' filesep 'RESULTS'];
cd(path);
files = dir([path '/*step*']);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
names={subFolders(:).name};
N=dir('*step*');
names={N(:).name};
names(contains(names(:, :), 'initial')) = [];
% names={'lungef1step1_I','lungef1step1_A','lungef1step1_P','lungef1step1_E'}
time=linspace(0, 1, 41);
Results=cell([length(names), 4]);
disp_rot=struct;

path=pwd;
for k = 1 : length(names)
    trialname=names{k};%(1:end-4);
    % trialname='lungef1step1_I';
    %file=[path filesep 'results.xlxs'];
    
    
     fig=  figure('visible', 'off');
     hold on
    
    if contains(trialname, '_A')==1 || contains(trialname, '_P')==1
        file=[path filesep trialname filesep 'target_AP.inp'];
        [fid,err] = fopen(file);
        if fid<=0 %does not exist
            disp(['file ' trialname ' does not exist'])
            continue
        else
            
            
            if all(fgetl(fid) == -1)
                % file is empty
                continue
            else
                l=dlmread(file);
                plot(time(1:length(l(:,1))), l(:,1));
                title(trialname)
                xlabel('time')
                ylabel('AP trans')
                displ=(l(end,1)-l(1,1));
                Results(k,1:2)= {trialname,num2str(displ)};
                disp_rot.(trialname)= num2str(displ);
            end
        end
    elseif  contains(trialname, '_I')==1 || contains(trialname, '_E')==1
        file=[path filesep trialname filesep 'target_IE.inp'];
          [fid,err] = fopen(file);
        if fid<=0 %does not exist
            continue
        else
            
            
            if all(fgetl(fid) == -1)
                % file is empty
                continue
            else
                l=dlmread(file);
                plot(time(1:length(l(:,1))), l(:,1));
                title(trialname)
                xlabel('time')
                ylabel('IE rot')
                rot=57.2957795*(l(end,1)-l(1,1));
                Results(k,3:4)= {trialname,num2str(rot)};
                  disp_rot.(trialname)= num2str(rot);
            end
        end
    end
     saveas(fig,[path filesep 'images' filesep trialname ],'jpeg')
     Results(~cellfun('isempty',Results)) ;
    
    
close('all');
    
end
   
%     if contains(trialname, '_A')==1 || contains(trialname, '_P')==1
%     file=[path filesep trialname filesep 'diag-ap.txt'];
%     else
%         file= [path filesep trialname filesep 'diag-ie.txt'];
%     end
%      fid = fopen(file);
%          if all(fgetl(fid) == -1)
%                 % file is empty
%                 continue
%          else
%                    data=dlmread(file);
%                    T=data(:,1);
%                    target=data(:,2); %-data(1,2);
%                    current=data(:,3);
%                    error=data(:,4);
%                    propErr=data(:,5);
%                    intErr=data(:,6);
%                     fig2=  figure('visible', 'off');
%              title(trialname)
%                 xlabel('time')
%                 subplot(2,1,1),plot(T, target, T, current);
%                   legend('Target', 'Current')
%                   subplot(2,1,2),
%                   yyaxis left
%                   ylabel('Proportional error')
%                   plot(T, propErr)
%                   yyaxis right
%                 plot(T,intErr);
%                     ylabel('Intergral error')
%                     legend('Prop error', 'Int error')
%                      saveas(fig2,[path filesep 'images' filesep 'diag' filesep trialname  '_diag'],'jpeg')
%                 close('all');
%          end
%         
% end


save([path  '\DISP_ROT.mat'],'disp_rot');
xlswrite([path filesep  'resultsall.xlsx'],Results);
fclose('all'); % close file


