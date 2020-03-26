close all

original_dir=pwd;
load('results.mat','RESULTS')
info_path= [ './jw/LtModel/force_amp_profiles/'];
load([info_path    '/Kneeflex.mat']);


impth=[pwd filesep 'jw' filesep 'RESULTS' filesep 'MUSCLEPLOTS'  ];
musc.A={  'gaslat' 'gasmed' ...
    'bflh' 'bfsh' 'grac'   'semimem' 'semiten' 'tfl' 'sart'};
%9
musc.P=    {    'recfem'  'vasint' 'vaslat' 'vasmed'};

musc.I=  {  'gaslat' ...
    'recfem'  'vasint' 'vaslat' 'vasmed' ...
    'bflh' 'bfsh'    'tfl' }; %External rotators lateral hammy

musc.E={   'gasmed' ...
    'recfem'  'vasint' 'vaslat' 'vasmed' ...
    'grac'   'semimem' 'semiten'  'sart'}; %internal rotators =medial

colours.E=createcolors(length( musc.E));

colours.A=createcolors(length( musc.A));

colours.I=createcolors(length( musc.I));

colours.P=createcolors(length( musc.P));

lines={'-', '--',':','-.'};

forcedata=[RESULTS.FORCE_MUSC];
muscleforces=struct;
% [o,inf]=structval(forcedata,'implant1')
% struct2cell(struct2cell(forcedata)(:))
% struct2table(forcedata)
% names = extractfield(forcedata,'implant1')
%
% types={'gait', 'lunge', 'staird', 'stairup', 'stepup'};
% nos=[1 2 3 4 5 ];
% step=1:4
%
% original_dir=pwd;

trials=fieldnames(forcedata);

lds=fieldnames(forcedata.lungef1.step1);

if isfield(forcedata.('lungef1'),'step2')==1
    forcedata.('lungef1')=rmfield(forcedata.('lungef1'),'step2');
    
end
if isfield(forcedata.ngait_og1.step3,'A')==1
        forcedata.ngait_og1.step3=rmfield(forcedata.ngait_og1.step3,'A');
        
    end

fs = 40; % data sample frequency
fc = 5;                        % cutoff frequency (arbitrary)
N = 4;                          % nth-order lowpass digital filter
Wn = fc/fs;

[b,a] = butter(N, Wn); %Butterworth IIR filter design
[bkin, akin] = butter(N, 0.7);


%muscles={'RF', 'VI','VM', 'VL','SMEMB', 'STEND', 'BFSH','BFLH'};
i=1; %count for anova
Anova=struct;
anova=struct;
ac=0;
c=1;
for imp=1:5
    impno=num2str(imp);
    %     i=0;
    APc=0; %counter without dividing into trials
    IEc=0;
    for idx1=1:length(trials)
        trial=trials{idx1};
        %trial count
        ac=ac+1;
        
        steps=fieldnames(forcedata.(trial));
        for idx2=1:length(steps)
            step=steps{idx2};
            sc=1; %step count
            
            lds=fieldnames(forcedata.(trial).(step));
            for idx3=1:length(lds)
                ld=lds{idx3};
                
                
                muscles=musc.(ld);
               % try
                    time=forcedata.(trial).(step).(ld).(['implant' impno])(2:end,1);
                    Tforce=zeros(1,length(time))';
                    for m=1:length(muscles)
                        muscle=muscles{m};
                        
                        musf=forcedata.(trial).(step).(ld).(['implant' impno])(2:end,2*m+1);
                        muscleforces.(trial).(step).(ld).(['imp' impno]).(muscle)=musf; %filtfilt(b,a,musf);
                        muscleforces.(trial).(step).(ld).time=time;
                        Tforce(:,1)=Tforce(:,1)+muscleforces.(trial).(step).(ld).(['imp' impno]).(muscle);
                        
                    end
                    muscleforces.(trial).(step).(ld).(['imp' impno]).total=Tforce;
                    muscleforces.totalall.data(c,1)=Tforce(end);
                    muscleforces.totalall.normdata(c,1)=Tforce(end)/totalforce.(trial).(step).(ld);
                    muscleforces.totalall.normdata2(c,1)=(Tforce(end)-totalforce.(trial).(step).(ld))/totalforce.(trial).(step).(ld);

                    %groups
                    muscleforces.totalall.implant(c,1)={['imp' impno]};
                    muscleforces.totalall.trial(c,1)={trial};
                    muscleforces.totalall.step(c,1)={step};
                    muscleforces.totalall.load(c,1)={ld};
                    c=c+1;
                    
%                     anova.(step).(ld).totalmuscleforce(:,ac)=Tforce(:,1);
%                     anova.(step).(ld).g1(ac)={trial};
%                     anova.(step).(ld).g2(ac)={['imp' impno]};
             %   catch
                  %  disp([ 'error with: ' trial ' ' step ' load' ld ' imp' impno])
              %  end
            end
        end
    end
end

Trials=(muscleforces.totalall.trial);
T2=struct2table(muscleforces.totalall)
TrialsCat=categorical(muscleforces.totalall.trial);

Implant=(muscleforces.totalall.implant);
ImplantCat=categorical(muscleforces.totalall.implant);

StepCat=categorical(muscleforces.totalall.step);
LoadCat=categorical(muscleforces.totalall.load);

D=muscleforces.totalall.data;
N=muscleforces.totalall.normdata;
N2=muscleforces.totalall.normdata2;

% T = table(D,N,Trials,Implant,StepCat,LoadCat);
Tstruct = table(D,N,TrialsCat,ImplantCat,StepCat,LoadCat);

AP=(Tstruct.LoadCat=='A' | Tstruct.LoadCat=='P');
TstructAP=Tstruct(AP,:);

IE=(Tstruct.LoadCat=='I' | Tstruct.LoadCat=='E');
TstructIE=Tstruct(IE,:);

I=(Tstruct.LoadCat=='I') ;
TstructI=Tstruct(I,:);

E=(Tstruct.LoadCat=='E') ;
TstructE=Tstruct(E,:);


%%
%Bar graph -trial type and implant. 


g2={'imp1','imp2','imp3','imp4','imp5'};
g1=trials';
barplot(Tstruct, g1, g2, 'Trials', 'Implant', 'Total muscle force comparing Implant and Trial type' )

barplot(TstructAP, g1, g2, 'Trials', 'Implant', 'Total muscle force for AP tests comparing Implant and Trial type')
barplot(TstructIE, g1, g2, 'Trials', 'Implant', 'Total muscle force for IE tests comparing Implant and Trial type')

g2={'imp1','imp2','imp3','imp4','imp5'};
g1=lds';
barplot(Tstruct, g1, g2, 'Load', 'Implant', 'Total muscle force comparing Implant and Load direction')

g2={'imp1','imp2','imp3','imp4','imp5'};
g1=steps';
barplot(Tstruct, g1, g2, 'Step', 'Implant', 'Total muscle force comparing Implant and Key instance')


% linearmodel = fitlm(Tstruct,'ResponseVar','D','PredictorVars',...
% {'TrialsCat', 'ImplantCat', 'LoadCat', 'StepCat'},'CategoricalVar',{'TrialsCat', 'ImplantCat', 'LoadCat', 'StepCat'})
% 
% linearmodel = fitrm(Tstruct,'ResponseVar','D','PredictorVars',...
% {'TrialsCat', 'ImplantCat', 'LoadCat'},'CategoricalVar',{'TrialsCat', 'ImplantCat', 'LoadCat'})
% 
% %%Anova
% meanByRegion = varfun(omean,T,'GroupingVariables','Region',...
%                       'InputVariables',{'Loss','Customers'})

[p,tab,stats] = anovan(table2array([Tstruct(:,2)]),{char(Tstruct.ImplantCat), char(Tstruct.LoadCat)...
      char(Tstruct.TrialsCat), char(Tstruct.StepCat) },...
 'varnames',{'Implant', 'Load', 'Trial','Step'} ) 


% 
% impl=table2cell(Td(:,4));', 'ImplantCat', 'LoadCat', 'StepCat'},...
% 'random',3,'varnames', {'Targets','HR','Subjects'});

   %         anova.(step).(ld).totalmuscleforce(end,:);
%         [p,~,stats] = anovan(datN,{tri impl},'varnames',{'Trial','Implant'},'display','off');
% 

% for idx2=1:length(steps)
%     step=steps{idx2};
%     for idx3=1:length(lds)
%         ld=lds{idx3};
%         
%         stepi=(T.StepCat)==step;
%         Td=T(stepi,:);
%         
%         ldi=((Td.LoadCat)==ld);
%         Td=Td(ldi,:);
%         dat=table2array(Td(:,1));
%         datN=table2array(Td(:,2));
%         tri=table2cell(Td(:,3));
%         impl=table2cell(Td(:,4));
 

  %% ANOVA 
      

% steps=    {'step1','step2','step3','step4'};
% 
% %compare between trial and implant
% for idx2=1:length(steps)
%     step=steps{idx2};
%     for idx3=1:length(lds)
%         ld=lds{idx3};
%         
%         stepi=(T.StepCat)==step;
%         Td=T(stepi,:);
%         
%         ldi=((Td.LoadCat)==ld);
%         Td=Td(ldi,:);
%         dat=table2array(Td(:,1));
%         datN=table2array(Td(:,2));
%         tri=table2cell(Td(:,3));
%         impl=table2cell(Td(:,4));
%         %         anova.(step).(ld).totalmuscleforce(end,:);
%         [p,~,stats] = anovan(datN,{tri impl},'varnames',{'Trial','Implant'},'display','off');
% 
%         statistics.(step).(ld).stat=stats;
%         pstat=cell(2,2);
%         pstat(1:2,1)=[{'Trial';'Implant'}];
%         pstat(1:2,2) =num2cell(p);
%         statistics.(step).(ld).p=pstat;
%         
%         [c,m,h]=multcompare(stats,'Dimension',[2 ]);
%         h.Name=['Multiple comparison of Trial type marginal means for ' step ' load ' ld];
%         saveas(h,[impth  filesep  'multcompareImplant_' step 'load' ld   ],'jpg');
%         
%         [c,m,h]=multcompare(stats,'Dimension',[1 ]);
%         h.Name=['Multiple comparison of implant type marginal means for ' step ' load ' ld];
%         saveas(h,[impth  filesep  'multcompareTrial_' step 'load' ld  ],'jpg');
%         
%         statistics.(step).(ld).stat=stats;
%     end
% end
% 
% % close all hidden

%%
%plot each muscle for each trial/step/load condition
impcolors={'b', 'y', 'r', 'g', 'm'};

for idx1=1:length(trials)
    trial=trials{idx1};
    if strcmp(trial,'lungef1')==1
        steps=    {'step1','step3','step4'};
    else
        steps=    {'step1','step2','step3','step4'};
    end
    for idx2=1:length(steps)
        step=steps{idx2};
        for idx3=1:length(lds)
            ld=lds{idx3};
            muscles=musc.(ld);
            fig= figure('visible', 'off');


% PatchInLegend = findobj(h_legend, 'type', 'patch');
            hold on
            for imp=1:5
                impno=num2str(imp);
                %                     try
                curves=[];
                for m=1:length(muscles)
                    muscle=muscles{m};
                    
                    
                    force=muscleforces.(trial).(step).(ld).(['imp' impno]).(muscle);
                    normfract=totalforce.(trial).fract.(step).(ld).(muscle);
                    
                    normF=(totalforce.(trial).(step).(ld));
                    Nforce=muscleforces.(trial).(step).(ld).(['imp' impno]).(muscle)/(normF*normfract);
                    curves=[force curves];
                    t= time;
                    % plot(t,force, 'Color',impcolors{imp}, 'LineWidth',1,'LineStyle',':')%colours.(ld)(m,:))%,'LineStyle',lines{imp})
                    
                    %plot individual muscles
                    
                    %                         plot(t,force, 'Color',impcolors{imp}, 'LineWidth',1,'LineStyle',':')%colours.(ld)(m,:))%,'LineStyle',lines{imp})
                    %                         title([trial ' ' step ' Load:' ld])
                    %                         xlabel('time');
                    %                         ylabel('Force (N)');
                end
 
                min_data=min(curves,[],2);
                max_data=max(curves,[],2);
                try
                plot(t,max_data, 'Color',impcolors{imp},'LineStyle','-', 'LineWidth',2')
                plot(t,min_data, 'Color',impcolors{imp},'LineStyle','-', 'LineWidth',2')
              
                h=patch([t;flipud(t)],[min_data;flipud(max_data)],impcolors{imp});
                set(h,'facealpha',.25)     ;
                catch
                    disp(['imp' impno ' ' trial ' ' step ' Load' ld 'failed'])
                end
                
            end
            L1 = plot(NaN,NaN,'b',NaN,NaN,'y',NaN,NaN,'r',NaN,NaN,'g',NaN,NaN,'m','LineWidth',2);
            [leg] = legend(L1, '0.5 imp1','0.4 imp2',  '0.3 imp3','0.2 imp4','0.1 imp5','Location','northwest');
            title(leg,'Conformity ratio')
            leg.Title.Visible = 'on';
            title([trial ' ' step ' Load:' ld])
            xlabel('time');
            ylabel('Force (N) ');
            
            
        
        
        saveas(fig,[impth  filesep  trial ' ' step ' load' ld  ],'jpg')
        close(fig)
        end
    end
end

function  barplot(T, group1, group2, group1label, group2label, Title)
Tstruct=T;
Bardata=struct;
tc=1;

Mean=[];
STD=[];

for idx1=1:length(group1)
    g1=group1{idx1};
    
    for idx2=1:length(group2)
        g2=group2{idx2};
        
        rows = Tstruct.([group1label 'Cat'])==g1 & Tstruct.([group2label 'Cat'])==g2;
        datatemp=table2array(Tstruct(rows,1));
        S = std(datatemp);
        M=mean(datatemp);
        Bardata.mean(tc)=M;
        Bardata.STD(tc)=S;
        Bardata.Load(tc)={g1};
        Bardata.Trial(tc)={g2};
        tc=tc+1;
        
        Mean(idx1,idx2)=M;
        STD(idx1,idx2)=S;
    end
end


figure
%g2=implant
ax = axes;
h = bar(Mean,'BarWidth',1);
h(1).FaceColor = 'b';
h(2).FaceColor = 'y';
h(3).FaceColor = 'r';
h(4).FaceColor = 'g';
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,1:length(group1));
% Naming each of the bar groups
xticklabels(ax,group1);
xlabel (group1label);
ylabel ('Force (N)');
title(Title)
% Creating a legend and placing it outside the bar plot
lg = legend(group2,'AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
ngroups = size(Mean, 1);
nbars = size(Mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Mean(:,i),STD(:,i), 'k', 'linestyle', 'none');
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

