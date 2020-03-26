for imp=[1:6]
implantno=num2str(imp);
cf = pwd;
pth=[cf filesep 'jw' filesep 'RESULTS_IMPLANT' implantno ];
files = dir([pth ]);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
n=size(subFolders,1);
nos=str2double({subFolders(:).name});
nos(isnan(nos)) = [];
nos=sort(nos);
folderno=nos(end);

PItimeplot({'all'},folderno,implantno)
end 




function [PIdata] = PItimeplot(trials,folderno,implantno, varargin)
if nargin==4 || nargin==5
    stepno=varargin{1};
    if nargin==5
        type=varargin{2};
    end
end
implantno=num2str(implantno);
currentFolder = pwd;
pth=[currentFolder filesep 'jw' filesep 'RESULTS_IMPLANT' implantno filesep num2str(folderno)];
impth=[currentFolder filesep 'jw' filesep 'IMAGES_IMPLANT' implantno ];

%addpath(genpath(pwd));
files = dir([pth '/*step*']);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
trialnames={subFolders(:).name};
trialnames(contains(trialnames(:, :), 'initial')) = [];

colors= createcolors(folderno+1);

if strcmp(trials{1},'all') ==1
else

    trialnames(~contains(trialnames(:, :), trials(:,:)))=[] ;
end
    if nargin==4 || nargin==5
        strn=cell(1,length(stepno)) ;
        for l=1:length(stepno)
            nos=strsplit(num2str(stepno));
            strn{l}=['step' nos{l}];
        end
        trialnames(~contains(trialnames(:, :),strn(:,:)))=[];
        
    end
    if nargin==5
        trialnames(~contains(trialnames(:, :), type)) = [];
    end

Results=cell([length(trialnames), 4]);

for i=1:length(trialnames)
    trialname=trialnames{i};
    % parentFolder=fileparts(pwd);
  
%     for k=1:folderno
%         pthN=[currentFolder filesep 'jw' filesep 'RESULTS_IMPLANT' implantno '_' num2str(k)];
%         if contains(trialname, '_A')==1 || contains(trialname, '_P')==1
%             fileN=[pthN filesep trialname filesep 'diag-ap.txt'];
%         else
%             fileN= [pthN filesep trialname filesep 'diag-ie.txt'];
%         end
%         [fidN,errN] = fopen(fileN, 'r');
%         if fidN<=0 %does not exist
%         
%             continue
%         else
%             dataN.(['trial' num2str(k)])=importdata(fileN);
%         end
%     end
    
    time=linspace(0, 1, 41);
    if contains(trialname, '_A')==1 || contains(trialname, '_P')==1
        file=[pth filesep trialname filesep 'target_AP.inp'];
        file2=[pth filesep trialname filesep 'diag-ap.txt'];
        Ylab='AP translation (mm)';
    else
        file2= [pth filesep trialname filesep 'diag-ie.txt'];
        file=[pth filesep trialname filesep 'target_IE.inp'];
        Ylab='IE rotation (degrees)';
    end
    [fid,err] = fopen(file, 'r');
    [fid2,err2] = fopen(file2, 'r');
    if fid<=0 ||  fid2<=0 %does not exist
        disp([ trialname ' file is does not exist'])
        continue
    else
        %         [fid,err] = fopen(file);
        %           [fid2,err2] = fopen(file2);
        
      
            
            if fseek(fid2, 1, 'bof') == -1
            disp([ trialname '  file is empty'])
            continue
        else
            
            fig=  figure('visible', 'off');
            hold on
            
            if fseek(fid, 1, 'bof') == -1
                l=zeros(1,41);
            else
            l=dlmread(file);
            end
           
            if contains(trialname, '_A')==1 || contains(trialname, '_P')==1
                 displ=(l(end,1)-l(1,1));
                Results(i,1:2)= {trialname,num2str(displ)};
            else
                 l=l*57.2957795;
                    displ=(l(end,1)-l(1,1));
                Results(i,3:4)= {trialname,num2str(displ)};
            end
            end
%             n=linecount(fid2);
%             if n<2
%                 disp([name ' empty']
%                 break
%             end
            % Read in file
% allData = textscan(fid2, '%s', 'delimiter', '\n');
% % Make allData cells empty if numerical 
% numericalArray = cellfun(@(s) sscanf(s,'%f').' ,allData, 'un', 0);
% % Get Header
% header = allData(cellfun('isempty',numericalArray));
% % Get Data
% data = vertcat(numericalArray{:});
           [ data]=importdata(file2);
%             data = cell2mat(f);
%             Z = textscan(fid2, '%f%f%f%f%f%f%f%f%f\n')%, 'delimiter', '\t'); % import full file
%             z = [Z{1}]; % separate into different cells
            fclose('all');
            
            T=data(:,1);
            target=data(:,2); %-data(1,2);
            current=data(:,3);
            error=data(:,4);
            propErr=data(:,5);
            intErr=data(:,6);
            
            PIdata.target.(trialname)=target;
            PIdata.current.(trialname)=current;
            PIdata.dt=(T(2)-T(1));
%             
%             subplot(3,1,1),
            plot(time(1:length(l))', l', 'g--o', T, target, T, current);
            title(['IMP' implantno ' ' trialname ' ' num2str(folderno)], 'Interpreter', 'none')
            xlabel('time')
           legend('Target', 'CurrTarget', 'Current','Location','eastoutside')


%        subplot(2,1,1), plot( time(1:length(l))', l', 'g--',T, target, T, current);
            

%             title(['IMP' implantno ' ' trialname ' ' num2str(folderno)], 'Interpreter', 'none')
%             
% ylabel(Ylab)
%            legend('Target', 'Current','Location','northwest')
% 
%             subplot(2,1,2)
%             plot(T, data(:,end))
%             ylabel('Total Muscle Force (N)')
% xlabel('time (s)')

%            
%             ylabel(Ylab)
%             
%             subplot(3,1,2)
%             yyaxis left
%             plot(T, propErr)
%             ylabel('Proportional error')
%             
%             yyaxis right
%             plot(T,intErr);
%             ylabel('Intergral error')
%             
%             legend('Prop error', 'Int error','Location','eastoutside')
%             
%             subplot(3,1,3)
%             plot(T, data(:,end))
%             ylabel('force')
%             if exist([impth  ],'dir')==0
%                 mkdir([impth  ])
%             end
            saveas(fig,[impth  filesep  trialname  '_' num2str(folderno) ],'jpg')
             close('all');
            
%             figN=figure('visible', 'off');
%             hold on
%             plot(time(1:length(l))', l', 'Color',colors(1,:), 'LineWidth', 2);
%             title([trialname ' all' ], 'Interpreter', 'none')
%             xlabel('time')
%             ylabel(Ylab)
%       
%             legendInfo{1} = 'target';
%       
%             for k=1:folderno
%                 legendInfo{k+1} = ['trial ' num2str(k)];
%                 if isfield(dataN,(['trial' num2str(k)]))==0
%                     continue
%                 else
%                     if isempty(dataN.(['trial' num2str(k)]))==1
%                         continue
%                     else
%                         T=dataN.(['trial' num2str(k)])(:,1);
%                         
%                         current=dataN.(['trial' num2str(k)])(:,3);
%                         
%                         
%                         plot(T, current,'Color',colors(k+1,:), 'LineWidth', 1);
%                     end
%                 end
%             end
%   
%            legend(legendInfo)
% 
%             saveas(figN,[impth  filesep  trialname  '_all' ],'jpg')
            close('all');
        end
    end
% end
Results(~cellfun('isempty',Results)) ;
%xlswrite([pth filesep  'resultsall.xlsx'],Results);
end


function n = linecount(fid)
n = 0;
tline = fgetl(fid);
while ischar(tline)
  tline = fgetl(fid);
  n = n+1;
end
end

function  colors= createcolors(n_colors)
 bg = [1 1 1];
 n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
    
    mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
    colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
end
    



