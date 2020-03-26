
function [finalloads] = totalmuscleloads
% Results=struct;
colors= createcolors(4);
% muscload=struct;
trialnames={'lungef1step1_A','lungef1step1_E','lungef1step1_I','lungef1step1_P',...
        'lungef1step3_A','lungef1step3_E','lungef1step3_I','lungef1step3_P',...
    'lungef1step4_A','lungef1step4_E','lungef1step4_I','lungef1step4_P',...
    'ngait_og1step1_A','ngait_og1step1_E','ngait_og1step1_I','ngait_og1step1_P',...
    'ngait_og1step2_A','ngait_og1step2_E','ngait_og1step2_I','ngait_og1step2_P',...
    'ngait_og1step3_A','ngait_og1step3_E','ngait_og1step3_I','ngait_og1step3_P',...
    'ngait_og1step4_A','ngait_og1step4_E','ngait_og1step4_I','ngait_og1step4_P',...
    'stairdownf2step1_A','stairdownf2step1_E','stairdownf2step1_I','stairdownf2step1_P',...
    'stairdownf2step2_A','stairdownf2step2_E','stairdownf2step2_I','stairdownf2step2_P',...
    'stairdownf2step3_A','stairdownf2step3_E','stairdownf2step3_I','stairdownf2step3_P',...
    'stairdownf2step4_A','stairdownf2step4_E','stairdownf2step4_I','stairdownf2step4_P',...
    'stairupf1step1_A','stairupf1step1_E','stairupf1step1_I','stairupf1step1_P',...
    'stairupf1step2_A','stairupf1step2_E','stairupf1step2_I','stairupf1step2_P',...
    'stairupf1step3_A','stairupf1step3_E','stairupf1step3_I','stairupf1step3_P',...
    'stairupf1step4_A','stairupf1step4_E','stairupf1step4_I','stairupf1step4_P',...
    'stepupf2step1_A','stepupf2step1_E','stepupf2step1_I','stepupf2step1_P',...
    'stepupf2step2_A','stepupf2step2_E','stepupf2step2_I','stepupf2step2_P',...
    'stepupf2step3_A','stepupf2step3_E','stepupf2step3_I','stepupf2step3_P',...
    'stepupf2step4_A','stepupf2step4_E','stepupf2step4_I','stepupf2step4_P'};
cf=pwd;
impth=[cf filesep 'jw' filesep 'IMAGES_MUSC' ];

fs = 40;
    % data sample frequency
    fc = 5;                        % cutoff frequency (arbitrary)
    N = 4;                          % nth-order lowpass digital filter
    Wn = fc/fs;
    
    [b,a] = butter(N, Wn); %Butterworth IIR filter design
    [bkin, akin] = butter(N, 0.7);
for imp=5
implantno=num2str(imp);
cf = pwd;
pth=[cf filesep 'jw' filesep 'RESULTS_IMPLANT' implantno ];
files = dir([pth ]);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
n=size(subFolders,1);
folderno=subFolders(n).name;

pthN=[cf filesep 'jw' filesep 'RESULTS_IMPLANT' implantno filesep folderno];


addpath(genpath(pwd));
files = dir([pthN '/*step*']);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
trialnames={subFolders(:).name};
trialnames(contains(trialnames(:, :), 'initial')) = [];


for i=1:length(trialnames)

    trialname=trialnames{i};
    % parentFolder=fileparts(pwd);

        if contains(trialname, '_A')==1 || contains(trialname, '_P')==1
            fileN=[pthN filesep trialname filesep 'diag-ap.txt'];
        else
            fileN= [pthN filesep trialname filesep 'diag-ie.txt'];
        end

        if exist(fileN,'file')==0 %does not exist
disp([trialname ' does not exisit'])
            continue
        else


            dataN=importdata(fileN);
              if isempty(dataN)==1
                  disp([trialname ' is empty'])
                        continue
                    else
            muscload.([trialname]).([ 'IMP' implantno])=filtfilt(b,a,dataN(:,end));
                        displace.([trialname]).([ 'IMP' implantno])=filtfilt(b,a,dataN(:,3));

              end
            

        end

end
 end
% % Results=cell([length(trialnames), 5]);
% 
 finalloads(:,1)=trialnames(:);
 save([cf filesep 'jw' filesep 'muscleloads.mat'],'muscload');
 save([cf filesep 'jw' filesep 'displace.mat'],'displace');%'-struct',


load([cf filesep 'jw' filesep 'muscleloads.mat'])
load([cf filesep 'jw' filesep 'displace.mat'])


for i=1:length(trialnames)
    fig=figure('visible','off');
    hold on
    trialname=trialnames{i};
    legendInfo={};
    title([trialname ' Muscle loads'], 'Interpreter', 'none')
    xlabel('time')
    for k=1:4
        implantno=num2str(k);
        
        if isfield(muscload.([trialname]),([ 'IMP' implantno]))==0
            legendInfo{k}=([' ']);
            continue
        else
            if isempty(muscload.([trialname]).([ 'IMP' implantno]))==1
                     legendInfo{k}=([' ']);
                continue
            else
                
                legendInfo{k} = (['IMP ' implantno]);
                current=muscload.([trialname]).([ 'IMP' implantno]);
                finalloads{i,k+1}=(current(end));
                T=linspace(0,1,size(current,1));
                switch k
                    case 1
                        c='r';
                    case 2
                        c='b';
                    case 3
                        c='y';
                    case 4
                        c='g';
                end
                plot(T, current,'Color',c, 'LineWidth', 1);
              
                
            end
        end
    end
    
    %legend(legendInfo)
       
    saveas(fig,[impth  filesep  trialname  ' muscleloads' ],'jpg')
    close('all');
end

for i=1:length(trialnames)
     trialname=trialnames{i};
     legendInfo={};
    figA=figure('visible','off');
    xlabel('time')
    title('Anterior')
     figP=figure('visible','off');
     xlabel('time')
      title('Posterior')
      figI=figure('visible','off');
      xlabel('time')
      title('Internal')
       figE=figure('visible','off');
       xlabel('time')
       title('External')
       
      y=displace.([trialname]).([ 'IMP5' ]);
      x=linspace(0,1,size(y,1));
       
       a=1;
       p=1;
       I=1;
       e=1;
       
       switch i
           case num2cell(1:12); color='r';
           case   num2cell(13:(12+16)); color='b';
                 case   num2cell((13+16):(12+2*16)); color='y';
                       case   num2cell((13+2*16):(12+3*16)); color='g';
                           
       end
    if contains(trialname, '_A')==1;        fig=figA;
         legendInfoA{a} = ([trialname]);
         set(0, 'CurrentFigure', figA)
         plot(x,y,'Color',color, 'LineWidth', 1);
     elseif  contains(trialname, '_P') ==1   ;      fig=figP;
         legendInfoP{p} = ([trialname]);
         set(0, 'CurrentFigure', figP)
         plot(x,y,'Color',color, 'LineWidth', 1);
    elseif  contains(trialname, '_I') ==1  ;      fig=figI;
        legendInfoI{I} = ([trialname]);
        set(0, 'CurrentFigure', figI)
        plot(x,y,'Color',color, 'LineWidth', 1);
    elseif  contains(trialname, '_E') ==1  ;      fig=figE;
        legendInfoE{e} = ([trialname]);
        set(0, 'CurrentFigure', figE)
        plot(x,y,'Color',color, 'LineWidth', 1);
        
                              
    end
     
    hold on
   
                 
                 
                
            end
    
    
    %legend(legendInfo)
       
    saveas(figA,[impth  filesep  'A_disp' ],'jpg')
    saveas(figP,[impth  filesep  'P_disp' ],'jpg')
    saveas(figI,[impth  filesep  'I_disp' ],'jpg')
    saveas(figE,[impth  filesep  'E_disp' ],'jpg')
    close('all');
end



% 
% Results(~cellfun('isempty',Results)) ;
% xlswrite([pth filesep  'resultsall.xlsx'],Results);




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

d={};

load('results.mat','RESULTS')
dispdata=[RESULTS.MOTIONS];

trials=fieldnames(forcedata);

lds=fieldnames(forcedata.lungef1.step1);

 dispdata.('lungef1')=rmfield(dispdata.('lungef1'),'step2');
lds=fieldnames(forcedata.lungef1.step1)
   steps=fieldnames(forcedata.stairupf1);
for q=1:4
    for imp=1:4
%     d{1,idx2*4-4+imp+1}=[step 'imp' num2str(imp)]      
        d{1,q*4-4+imp+1}=['step' num2str(q) 'imp' num2str(imp)]
    end
end

i=1;

for idx1=1:length(trials)
    trial=trials{idx1};
    for l=1:length( lds)
        ld=lds{l}
        
        switch ld
            case 'A'; c=2; n=1;
            case 'P'; c=2;n=1;
            case 'I'; c=3;n=57.2958;
            case 'E'; c=3;n=57.2958;
        end
        
        
        d{i+1,1}=[trial  ld]
        
        
        for idx2=1:length(steps)
            step=steps{idx2};
            for imp=1:4
                impno=num2str(imp);
                
                
                if isfield (dispdata.(trial),(step))==1
                    if isfield (dispdata.(trial).(step).(ld),(['implant' impno]))==1
                        d{i+1, idx2*4-4+imp+1}=n*dispdata.(trial).(step).(ld).(['implant' impno])(end,c)
                    else
                        continue;
                    end
                    else
                        continue
                    end
                end
            end
      
        i=i+1;
    end
end
