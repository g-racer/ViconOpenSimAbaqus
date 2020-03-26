close all hidden
original_dir=pwd;
load('results.mat','RESULTS')
info_path= [ './jw/LtModel/force_amp_profiles/'];
load([info_path    'kneeflex.mat']);
load([info_path    'GRF.mat']);

py={'FORCE_MUSC','CPRESS','CAREA','CFN','COP','JOINT_LOADS'};

vars.CPRESS={'TF_MED_PRESS' 'TF_LAT_PRESS'  'PF_PRESS'};
vars.COP={'MED_XN1' 'MED_XN2' 'MED_XN3'  'LAT_XN1' 'LAT_XN2' 'LAT_XN3'...
    'PAT_XN1' 'PAT_XN2' 'PAT_XN3'};
vars.CAREA={'MED_TF_AREA'     'LAT_TF_AREA'    'PF_AREA'};

vars.CPRESS_P={'TOTAL_CPRESS'  'PF_PRESS'};
vars.COP_P={'TF_XN1' 'TF_XN2' 'TF_XN3'  ...
    'PAT_XN1' 'PAT_XN2' 'PAT_XN3'};
vars.CAREA_P={'TOTAL_CAREA'       'PF_AREA'};
vars.JOINT_LOADS={'TF_CTF1' 'TF_CTM1'};
vars.JOINT_LOADS_P=vars.JOINT_LOADS;
vars.CFN={'PAT_CFN1'  'PAT_CFN2'  'PAT_CFN3'};
vars.CFN_P={'PAT_CFN1'  'PAT_CFN2'  'PAT_CFN3'};

impth=[pwd filesep 'jw' filesep 'RESULTS' filesep 'MUSCLEPLOTS'  ];
vars.FORCE_MUSC.A={  'gaslat' 'gasmed' ...
    'bflh' 'bfsh' 'grac'   'semimem' 'semiten' 'tfl' 'sart'};
%9
vars.FORCE_MUSC.P=    {    'recfem'  'vasint' 'vaslat' 'vasmed'};

vars.FORCE_MUSC.I=  {  'gaslat' ...
    'recfem'  'vasint' 'vaslat' 'vasmed' ...
    'bflh' 'bfsh'    'tfl' }; %External rotators lateral hammy

vars.FORCE_MUSC.E={   'gasmed' ...
    'recfem'  'vasint' 'vaslat' 'vasmed' ...
    'grac'   'semimem' 'semiten'  'sart'}; %internal rotators =medial

lines={'-', '--',':','-.'};
RESULTSstruct=struct;

figures = [];
% Generate figures
%[your code]
% Resize and output figures
figSize = [20, 20];            % [width, height]
figUnits = 'Centimeters';

for d=1:length(py)
    
    Data=[RESULTS.(py{d})];
    
    c=1;
    C=1;
    
    trials=fieldnames(Data);
    lds=fieldnames(Data.lungef1.step1);
    
    if isfield(Data.('lungef1'),'step2')==1
        Data.('lungef1')=rmfield(Data.('lungef1'),'step2');
    end
    if isfield(Data.('ngait_og1').step3,'A')==1
        Data.('ngait_og1').step3=rmfield(Data.('ngait_og1').step3,'A');
    end
    
    if isfield(Data.('ngait_og1').step3,'E')==1
        Data.('ngait_og1').step3=rmfield(Data.('ngait_og1').step3,'E');
    end
    
    conformity=[(5:-1:1)*0.1 0.6];
    for imp=[ 1:5 6]
        impno=num2str(imp);
        APc=0; %counter without dividing into trials
        IEc=0;
        for idx1=1:length(trials)
            trial=trials{idx1};
            steps=fieldnames(Data.(trial));
            for idx2=1:length(steps)
                step=steps{idx2};
                sc=1; %step count
                
                lds=fieldnames(Data.(trial).(step));
                for idx3=1:length(lds)
                    ld=lds{idx3};
                    
                    if contains(py{d},'FORC')==1
                        colheaders=vars.FORCE_MUSC.(ld);
                    else
                        if impno == '6' || impno=='7'
                            colheaders=vars.([py{d} '_P']);
                        else
                            colheaders=vars.(py{d});
                        end
                    end
                    % try
                    time=Data.(trial).(step).(ld).(['implant' impno])(2:end,1);
                    if contains(py{d},'FORC')==1
                        Tforce=zeros(1,length(time))';
                        % find max for east test to normalise
                        s=sum(RESULTS.FORCE_MUSC.(trial).(step).(ld).(['implant5'])(:,3:2:end),2);
                        Max=s(end);
                        
                        for m=1:length(colheaders)
                            header=colheaders{m};
                            
                            fill=Data.(trial).(step).(ld).(['implant' impno])(2:end,2*m+1);
                            RESULTSstruct.(py{d}).(trial).(step).(ld).(['imp' impno]).(header)=fill; %filtfilt(b,a,musf);
                            
                        end
                        Tforce=sum(Data.(trial).(step).(ld).(['implant' impno])(:,3:2:end),2);
                        
                        RESULTSstruct.(py{d}).(trial).(step).(ld).(['imp' impno]).total=Tforce;
                        RESULTSstruct.(py{d}).totalall.implant(c,1)={['imp' impno]};
                        RESULTSstruct.(py{d}).totalall.trial(c,1)={trial};
                        RESULTSstruct.(py{d}).totalall.step(c,1)={step};
                        RESULTSstruct.(py{d}).totalall.load(c,1)={ld};
                        RESULTSstruct.(py{d}).totalall.conf(c,1)=conformity(imp);
                        %                        RESULTSstruct.(py{d}).totalall.flex(c,1)=Kneeflex.(trial).(step);
                        RESULTSstruct.(py{d}).totalall.GRF(c,1)=GRF.(trial).(step);
                        RESULTSstruct.(py{d}).totalall.FORCE_MUSC(c,1)=Tforce(end);
                        
                        RESULTSstruct.(py{d}).totalall.FORCE_MUSC_Norm(c,1)=Tforce(end)/Max;% totalforce.(trial).(step).(ld);
                        %       RESULTSstruct.(py{d}).totalall.FORCE_MUSC_Nflex(c,1)=Tforce(end)/kneeflex.(trial).(step);
                        
                        c=c+1;
                    else
                        for m=1:length(colheaders)
                            header=colheaders{m};
                            fill=Data.(trial).(step).(ld).(['implant' impno])(2:end,m+1);
                            RESULTSstruct.(py{d}).(trial).(step).(ld).(['imp' impno]).(header)=fill; %filtfilt(b,a,musf);
                            
                            %overtime
                            RESULTSstruct.(py{d}).(trial).(step).(ld).(['imp' impno]).(header)=fill; %filtfilt(b,a,musf);
                            %end
                            RESULTSstruct.(py{d}).end.(header)(C,1)=Data.(trial).(step).(ld).(['implant' impno])(end,m+1);
                        end
                        if imp~=6 && imp~=7
                            if contains(py{d},'CPRESS')==1
                                RESULTSstruct.CPRESS.end.TOTAL_CPRESS(C,1)=RESULTSstruct.CPRESS.end.TF_MED_PRESS(C,1)+RESULTSstruct.CPRESS.end.TF_LAT_PRESS(C,1);
                            end
                            if contains(py{d},'CAREA')==1
                                RESULTSstruct.CAREA.end.TOTAL_CAREA(C,1)=RESULTSstruct.CAREA.end.MED_TF_AREA(C,1)+RESULTSstruct.CAREA.end.LAT_TF_AREA(C,1);
                            end
                            if contains(py{d},'COP')==1
                                RESULTSstruct.COP.end.TF_XN1(C,1)=(RESULTSstruct.COP.end.MED_XN1(C,1)+RESULTSstruct.COP.end.LAT_XN1(C,1))/2;
                                RESULTSstruct.COP.end.TF_XN2(C,1)=(RESULTSstruct.COP.end.MED_XN2(C,1)+RESULTSstruct.COP.end.LAT_XN2(C,1))/2;
                                RESULTSstruct.COP.end.TF_XN3(C,1)=(RESULTSstruct.COP.end.MED_XN3(C,1)+RESULTSstruct.COP.end.LAT_XN3(C,1))/2;
                                
                            end
                        end
                        
                        RESULTSstruct.(py{d}).(trial).(step).(ld).time=time;
                        %final values
                        RESULTSstruct.(py{d}).end.implant(C,1)={['imp' impno]};
                        RESULTSstruct.(py{d}).end.trial(C,1)={trial};
                        RESULTSstruct.(py{d}).end.step(C,1)={step};
                        RESULTSstruct.(py{d}).end.load(C,1)={ld};
                        
                        C=C+1;
                    end
                    
                end
                
            end
        end
        
    end
    
end

%% create tables

rn=fieldnames(RESULTSstruct.CPRESS.end);
fields=rn(contains(rn(:, :), 'TF'));% = [];
CP = rmfield(RESULTSstruct.CPRESS.end,fields);
CPRESS=struct2table(CP);


rn=fieldnames(RESULTSstruct.COP.end);
fields=rn(contains(rn(:, :), 'MED'));% = [];
CoP = rmfield(RESULTSstruct.COP.end,fields);
fields=rn(contains(rn(:, :), 'LAT'));% = [];
CoP = rmfield(CoP,fields);
COP=struct2table(CoP);

rn=fieldnames(RESULTSstruct.CAREA.end);
fields=rn(contains(rn(:, :), 'TF'));% = [];
CA = rmfield(RESULTSstruct.CAREA.end,fields);
CAREA=struct2table(CA);


CFN=struct2table(RESULTSstruct.CFN.end);
CFN.PAT_CFN=sqrt(CFN.PAT_CFN1.^2+CFN.PAT_CFN2.^2+CFN.PAT_CFN3.^2);

JOINT_LOADS=struct2table(RESULTSstruct.JOINT_LOADS.end);

MUSC_FORCE=struct2table(RESULTSstruct.FORCE_MUSC.totalall);
ResultsTable=struct2table(RESULTSstruct.FORCE_MUSC.totalall);
% for d=2:length(py)
%
%     ResultsTable = innerjoin(ResultsTable,struct2table(RESULTSstruct.(py{d}).end(:,1)),'Keys',{'implant' 'trial' 'step' 'load'});
% end

ResultsTable = innerjoin(ResultsTable,CAREA,'Keys',{'implant' 'trial' 'step' 'load'});
ResultsTable = innerjoin(ResultsTable,CPRESS,'Keys',{'implant' 'trial' 'step' 'load'});
ResultsTable = innerjoin(ResultsTable,JOINT_LOADS,'Keys',{'implant' 'trial' 'step' 'load'});
ResultsTable.TF_CTM1=abs(ResultsTable.TF_CTM1);
%  ResultsTable = innerjoin(ResultsTable,MUSC_FORCE,'Keys',{'implant' 'trial' 'step' 'load'});


ResultsStruct=table2struct(ResultsTable) ;
%
%%
% depvar={'FORCE_MUSC'  'FORCE_MUSC_Nflex' 'FORCE_MUSC_Norm' 'ML_PRESS'  'ML_AREA' 'MED_TF_AREA' 'LAT_TF_AREA' 'PF_AREA'...
%     'PAT_CFN1'  'PAT_CFN2'  'PAT_CFN3'};
depvar={'FORCE_MUSC'  'FORCE_MUSC_Norm' 'PF_AREA' 'TOTAL_CAREA' 'PF_PRESS' 'TOTAL_CPRESS' 'TF_CTF1' 'TF_CTM1' };
depvarlab={'Muscle Force'  'Normalised Muscle Force' 'PF Contact Area' 'Total Contact Area'...
    'PF Contact Pressure' 'Total Contact Pressure' 'Axial Joint Load' 'Axial Joint Moment' };

inpvar={'implant' 'trial' 'step' 'load'};
dsa=ResultsTable(:,{inpvar{:},depvar{:}});

ylabels={'Force (N)'  'Force/Max Force ' 'Area (mm^2)' 'Area (mm^2)'  'Pressure (Pa)' 'Pressure (Pa)' 'Force (N)' 'Moment (Nm)'};
columns=[7:14];


% ylabels={'FORCE (N)' 'FORCE/FORCE ' 'FORCE/FORCE ' 'PRESSURE (Pa)' 'PRESSURE (Pa)' 'PRESSURE (Pa)' ...
%     'AREA (mm^2)' 'AREA (mm^2)' 'AREA (mm^2)','FORCE (N)','FORCE (N)','FORCE (N)'};
id1=find(strcmp(ResultsTable.load, 'A'));
id2=find(strcmp(ResultsTable.load, 'P'));
id3=union(id1, id2);

ResultsTableAP=ResultsTable(id3,:);

id1=find(strcmp(ResultsTable.load, 'I'));
id2=find(strcmp(ResultsTable.load, 'E'));
id3=union(id1, id2);

ResultsTableIE=ResultsTable(id3,:);

%%

for id=1:length(depvar)
    ylab=ylabels{id} ;
    %     %statarray = grpstats(dsa,depvar{id},{'mean','std','gName'})
    %
    %

    %ALL tests
    [p,tab,stats] = anovan(ResultsTable.(depvar{id})(:),{char(ResultsTable.implant), char(ResultsTable.trial)},...
        'varNames',{'Implant',  'Trial'},'display','off');
    figure;
    [c,meadstd,handle,grp]=multcompare(stats, 'Dimension',[ 1 2],'CType','bonferroni');
    multcmp.(['stats' depvar{id}]).groupcomparisonimplant =c ;
    multcmp.(['stats' depvar{id}]).pvalue =[{'Implant';  'Trial'} num2cell(p)] ;
    %multcmp.(['statsAP' depvar{id}]).meanSTD =[{'Implant1 -0.5',  'Implant2 -0.4','Implant3 -0.3','Implant4 -0.2','Implant5 -0.1'}' num2cell(meadstd)] ;
    multcmp.(['stats' depvar{id}]).meanSTD =[grp num2cell(meadstd)] ;
    
    handle.Name=['Multiple comparison of population marginal means for ' (depvarlab{id}) ];
    
    % %     if h==0
    
%     %AP tests
%     [p,tab,statsAP] = anovan(ResultsTableAP.(depvar{id})(:),{char(ResultsTableAP.implant), char(ResultsTableAP.trial)},...
%         'varNames',{'Implant',  'Trial'},'display','off');
%     figure;
%     [cAP,meadstd,handle,grp]=multcompare(statsAP, 'Dimension',[ 1 2]);
%     multcmp.(['statsAP' depvar{id}]).groupcomparisonimplant =cAP ;
%     multcmp.(['statsAP' depvar{id}]).pvalue =[{'Implant';  'Trial'} num2cell(p)] ;
%     %multcmp.(['statsAP' depvar{id}]).meanSTD =[{'Implant1 -0.5',  'Implant2 -0.4','Implant3 -0.3','Implant4 -0.2','Implant5 -0.1'}' num2cell(meadstd)] ;
%     multcmp.(['statsAP' depvar{id}]).meanSTD =[grp num2cell(meadstd)] ;
%     
%     handle.Name=['Multiple comparison of population marginal means for ' (depvarlab{id}) ' AP test'];
%     
%     % anova.(['statsAP' depvar{id}]) =p;
%     
%     %IE tests
%     [p,tab,statsIE] = anovan(ResultsTableIE.(depvar{id})(:),{char(ResultsTableIE.implant), char(ResultsTableIE.trial)}...
%         , 'varNames',{'Implant',  'Trial'},'display','off') ;
%     anova.(['statsIE' depvar{id}]) =p;
%     figure;
%     [cIE,meadstd,handle,grp]=multcompare(statsIE, 'Dimension',[ 2 1]);
%     multcmp.(['statsIE' depvar{id}]).pvalue =[{'Implant';  'Trial'} num2cell(p)] ;
%     multcmp.(['statsIE' depvar{id}]).meanSTD =[grp num2cell(meadstd)] ;
%     multcmp.(['statsIE' depvar{id}]).Groupcomparison =cIE ;
%     handle.Name=['Multiple comparison of population marginal means for ' (depvarlab{id}) ' IE test'];
%     
    %       [p,table,stats] = anova1(ResultsTableAP.FORCE_MUSC_Nflex(:),(ResultsTableAP.implant(:) ));
    %
    %
    %     else
    % [p,tbl,stats] = kruskalwallis(ResultsTable.(depvar{id})(:),{char(ResultsTable.implant), char(ResultsTable.trial)...
    %         char(ResultsTable.step), char(ResultsTable.load) });%,...
    %  'varNames',{'Implant',  'Trial','Step','Load'} ) ;
    %     end
    
    
    g2={ 'imp1','imp2','imp3','imp4','imp5', 'imp6'};
    g1=categories(categorical(ResultsTable.trial))';
      [X, Y,E,P,grps, Pr, Prgrps] =  barsetup(ResultsTable(:,[1:4 columns(id)]), g1, g2, 'trial', 'implant', c, grp);
    groupsig.(depvar{id})=grps;
    Cl = nan(5, 5, 3);
 
    Cl(1, :, :) = repmat([.34 .63 .83],5,1);
    Cl(2, :, :) = repmat([.6 .4 .8],5,1);
    Cl(3, :, :) = repmat([.3 .9 .3],5,1);
    Cl(4, :, :) = repmat([0.88, 0.52, 0.35],5,1);
    Cl(5, :, :) = repmat([1 .94 .0],5,1);
    Cl(6, :, :) = repmat([.18 .52 .49],5,1);
    C=reshape(Cl(:,1,:),6,3);
    %
    Colors= [0, 0.4470, 0.7410;
        0.9290, 0.6940, 0.1250;
        0.8500, 0.3250, 0.0980;
        0.4660, 0.6740, 0.1880;
        0.4940, 0.1840, 0.5560  ];
    h0=figure;
    h0.Name=[depvarlab{id} ' Comparing Implant and Trial type'];
    hold on
    lab= {'Lunge' 'Gait','Stairdown','Stairup','Stepup'};
    %superbar(X',Y', 'E', E', 'P', P, 'BarFaceColor', Cl,'PStarShowNS',false);
    superbar(X',Y', 'P', Pr, 'E', E', 'BarFaceColor', Cl,...
        'PStarShowNS',false,'PStarFontSize',8,'PStarIcon',char('(L)'));
    tit=sprintf([depvarlab{id} ' Comparing\n Implant and Trial type']);
    figuremaker(lab,  'Trial',  tit , ylab,C)
    savefigure(h0);

    
    %%
%     [X, Y,E,P,grps, PrAP, Prgrps] =  barsetup(ResultsTableAP(:,[1:4 columns(id)]), g1, g2, 'trial', 'implant', cAP, grp);
%     groupsigAP.(depvar{id})=grps;
%    
%     h0=figure;
%     h0.Name=[depvarlab{id} ' Comparing Implant and Trial type -AP'];
%     hold on
%     lab= {'Lunge' 'Gait','Stairdown','Stairup','Stepup'};
%     %superbar(X',Y', 'E', E', 'P', P, 'BarFaceColor', Cl,'PStarShowNS',false);
%     superbar(X',Y', 'P', PrAP, 'E', E', 'BarFaceColor', Cl,...
%         'PStarShowNS',false,'PStarFontSize',8,'PStarIcon',char('(L)'));
%     tit=sprintf([depvarlab{id} ' Comparing\n Implant and Trial type -AP']);
%     figuremaker(lab,  'Trial',  tit , ylab,C)
%     savefigure(h0);
%     
%     
%     %     figures=[figures h0];
%     [X, Y,E,P, grps, PrIE, Prgrps] =  barsetup(ResultsTableIE(:,[1:4 columns(id)]), g1, g2, 'trial', 'implant', cIE, grp);
%     h0=figure;
%     groupsigIE.(depvar{id})=grps;
%     h0.Name=[depvarlab{id} ' Comparing Implant and Trial type -IE'];
%     hold on
%     % superbar(X',Y', 'E', E', 'P', P, 'BarFaceColor', Cl,'PStarShowNS',false);
%     superbar(X',Y', 'E', E', 'P', PrIE, 'BarFaceColor', Cl,...
%         'PStarShowNS',false,'PStarFontSize',8,'PStarIcon',char('(L)'));
%     tit=sprintf([depvarlab{id} ' Comparing\n Implant and Trial type -IE']);
%     figuremaker(lab,  'Trial',  tit , ylab,C)
%     %figures=[figures h0];
%     savefigure(h0);
    %%
    
    %     barplot(ResultsTableAP(:,[1:4 columns(id)]), g1, g2, 'trial', 'implant', [depvar{id} ' Comparing Implant and Trial type -AP'] , ylab,cAP);
    %     barplot(ResultsTableIE(:,[1:4 columns(id)]), g1, g2, 'trial', 'implant', [depvar{id} ' Comparing Implant and Trial type -IE'] , ylab, cIE);
    %     barplot(ResultsTable(:,[1:4 columns(id)]), g1, g2, 'trial', 'implant', [depvar{id} ' Comparing Implant and Trial type '] , ylab)
    %
    %
    
    %     barplot(TstructAP, g1, g2, 'Trials', 'Implant', 'Total muscle force for AP tests Comparing Implant and Trial type')
    %     barplot(TstructIE, g1, g2, 'Trials', 'Implant', 'Total muscle force for IE tests Comparing Implant and Trial type')
    %
    %     g2={'imp1','imp2','imp3','imp4'};
    %     g1=lds';
    %     barplot(Tstruct, g1, g2, 'Load', 'Implant', 'Total muscle force Comparing Implant and Load direction')
    %
    %     g2={'imp1','imp2','imp3','imp4'};
    %     g1=steps';
    %     barplot(Tstruct, g1, g2, 'Step', 'Implant', 'Total muscle force Comparing Implant and Key instance')
    %
end

g1=categories(categorical(ResultsTable.trial))';
barplot1(ResultsTable(:,[1:4 9]),  g2,  'implant', ['FORCE_MUSC_Norm  Comparing Implant and Trial type '] ,'FORCE/FORCE ');

% bar(categorical(ResultsTable.implant),ResultsTable.FORCE_MUSC_Norm) ['FORCE_MUSC_Nflex Comparing Implant and Trial type '] , ylab);

%% regression %%
%for all
Resultsregress.all=ResultsTable;
id6=find(strcmp(Resultsregress.all.implant, 'imp6'));
id7=find(strcmp(Resultsregress.all.implant, 'imp7'));
id=union(id6, id7);
Resultsregress.all(id,:)=[];


depvar={ 'FORCE_MUSC_Norm' 'PF_AREA' 'TOTAL_CAREA' 'PF_PRESS' 'TOTAL_CPRESS' 'TF_CTF1' };
depvarlab={  'Normalised Muscle Force' 'PF Contact Area' 'Total Contact Area'...
    'PF Contact Pressure' 'Total Contact Pressure' 'Axial Joint Load'  };
ylabels={ ['Force/Max Force '] ['Area (mm^{2})'] ['Area (mm^{2})']  ['Pressure (Pa)'] ['Pressure (Pa)'] ['Force (N)'] };

inpvar={'implant' 'trial' 'step' 'load'};
for id=1:length(depvar)
    ylab=ylabels{id} ;
    label=depvarlab{id};
    variable=depvar{id};

h2= figure;
h2.Name=['Plot of ' label ' for all Loads'];
% regress=Resultsregress.all;
%     X=regress.conf;
%     Y=regress.(variable);
%     modelspec ='linear';
%     mdl = fitlm(X,Y);
    

hold on
R2=cell(length(trials),1);
for k=1:length(trials)
    trial=trials{k};
    
    id1=find(strcmp(Resultsregress.all.trial, trial));
    
    regress=Resultsregress.all(id1,:);
    X=regress.conf;
    Y=regress.(variable);
    modelspec ='linear';
    mdl = fitlm(X,Y);
    %
%  m = 2; b = 1;
% fplot(@(x)m*x+b, [0.1 0.5]);
%         fig= plot(linspace(0.1,0.5,size(why,1)),why);

    fig= plot(mdl);
    fig(1).Color=imcomplement(Colors(k,:));%Colors(k,:);
    fig(2).Color=imcomplement(Colors(k,:));
    fig(3).Color='none';%[imcomplement(Colors(k,:))];%+Colors(k,:)]/2;
    fig(4).Color='none';%[imcomplement(Colors(k,:))];%+Colors(k,:)]/2;
    
    for g=1:4
        fig(g).LineWidth=2;
        fig(g).MarkerSize=8;
    end
    %temp=%num2str(mdl.Rsquared.Ordinary);
    R2{k,1}=num2str(sprintf('%1.3f\n',mdl.Rsquared.Ordinary));

 %   R2{k,1}=mdl.Rsquared.Ordinary;
end
%title(['Plot of Normalized Muscle Force for all Loads']);
xticks(0:0.1:0.5)
title('')
xlabel('Conformity Ratio');
ylabel([ylab],'interpreter','tex');
L1=plot(NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,'LineWidth',2);
for i=1:5
    L1(i).Color=imcomplement(Colors(i,:));
end

[leg] = legend(L1, ['Lunge R^{2}=' R2{1}], [ 'Gait R^{2}=' R2{2}] , ['Stairdown R^{2}=' R2{3}],...
    ['Stairup R^{2}=' R2{4}], ['Stepup R^{2}=' R2{5}],'Location','northeast');
title(leg,'Trial')
leg.Title.Visible = 'on';
savefigure(h2);
%figures=[figures h2];
%for each trial
end
%%
idA=find(strcmp(Resultsregress.all.load, 'A'));
idP=find(strcmp(Resultsregress.all.load, 'P'));
idAP=union(idA, idP);

idI=find(strcmp(Resultsregress.all.load, 'I'));
idE=find(strcmp(Resultsregress.all.load, 'E'));
idIE=union(idI, idE);

for id=1:length(depvar)
    ylab=ylabels{id} ;
    label=depvarlab{id};
    variable=depvar{id};
    
    for l=1:2
        if l==1
            ind=idAP;
            ld='AP';
        else
            ind=idIE;
            ld='IE';
        end
        
        %
        %     id2=find(strcmp(Resultsregress.all.load, ld));
        Resultsregress.(ld)=Resultsregress.all(ind,:);
        
        h3= figure;
        h3.Name=['Plot of ' label ' for ' ld ' Loads'];
        hold on
        
        for k=1:length(trials)
            trial=trials{k};
            
            id1=find(strcmp(Resultsregress.(ld).trial, trial));
            
            regress=Resultsregress.(ld)(id1,:);
            X=regress.conf;
            Y=regress.(variable);
            modelspec ='linear';
            mdl = fitlm(X,Y);
            %
            
            fig= plot(mdl);
            fig(1).Color=imcomplement(Colors(k,:));%Colors(k,:);
            fig(2).Color=imcomplement(Colors(k,:));
            fig(3).Color='none';%[imcomplement(Colors(k,:))];%+Colors(k,:)]/2;
            fig(4).Color='none';%[imcomplement(Colors(k,:))];%+Colors(k,:)]/2;
            for g=1:4
                fig(g).LineWidth=2;
                fig(g).MarkerSize=8;
            end
            
            R2{k,1}=num2str(sprintf('%1.3f\n',mdl.Rsquared.Ordinary));
            
        end
        %title(['Plot of Normalized Muscle Force for all Loads']);
        xticks(0:0.1:0.5)
        title('')
        xlabel('Conformity Ratio');
ylabel([ylab],'interpreter','tex');
        L1=plot(NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,'LineWidth',2);
        for i=1:5
            L1(i).Color=imcomplement(Colors(i,:));
        end
        
        [leg] = legend(L1, ['Lunge R^{2}=' R2{2}], [ 'Gait R^{2}=' R2{1}] , ['Stairdown R^{2}=' R2{3}],...
            ['Stairup R^{2}=' R2{4}], ['Stepup R^{2}=' R2{5}],'Location','northeast');
        title(leg,'Trial')
        leg.Title.Visible = 'on';
        
        savefigure(h3);
        
    end
end

    % figures=[figures h3];

%%
function savefigure(fig)
% for f = 1:numel(figures)
% fig = figures(f);
% Resize the figure
%       set(fig, 'Units', figUnits);
%       pos = get(fig, 'Position');
%       pos = [pos(1), pos(4)+figSize(2), pos(3)+figSize(1), pos(4)];
%       set(fig, 'Position', pos);
Name=fig.Name;
fileName = sprintf('%s', [Name '.jpeg']);

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
% print(fig, '-dpdf', 'test3.pdf', '-append');
% Output the figure
%  print(fig, '-dpsc', '-append', fileName)
print( fig, '-djpeg', fileName);
% end
end
% %regression analysis
% modelspec = 'FORCE_MUSC ~ conf*GRF*flex - conf:GRF:flex';
% mdl = fitlm(ResultsTable,modelspec)
% y=ResultsTable.FORCE_MUSC;
% ResultsTable.C=categorical(ResultsTable.conf);
%
% %
% x1=ResultsTable.flex;
% x2=ResultsTable.GRF;
% x3=ResultsTable.conf;
% X=[1./x1 ]%1./x2]% x3];
% mdl = fitlm(X,y)
% mdl = polyfit(x1,y,n)
% mdl = fitlm(ResultsTable,'FORCE_MUSC ~ conf + flex^2 + conf*flex' )%+ GRF + flex')
% mdl = fitnlm(ResultsTable,'FORCE_MUSC ~ conf ', 1 )%
% plot(mdl)



% for id=1:length(inpvar)
%     dsa=ResultsTable(:,{inpvar{id},depvar{:}});
%     statarray = grpstats(dsa,inpvar{id},{'mean','std','gName'});
% end
%





function  barplot1(T, group1, group1label, Title, ylab)
Tstruct=T;
Bardata=struct;
tc=1;

Mean=[];
STD=[];

for idx1=1:length(group1)
    g1=group1{idx1};
    
    rows = strcmp(Tstruct.([group1label]),g1)==1 ;
    datatemp=table2array(Tstruct(rows,5));
    S = std(datatemp);
    M=mean(datatemp);
    Bardata.mean(tc)=M;
    Bardata.STD(tc)=S;
    %         Bardata.Load(tc)={g1};
    %         Bardata.Trial(tc)={g2};
    tc=tc+1;
    
    Mean(idx1)=M;
    STD(idx1)=S;
    
end

figure
%g2=implant
ax = axes;
h = bar(Mean,'BarWidth',1);
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,1:length(group1));
% Naming each of the bar groups
xticklabels(ax,group1);
xlabel (group1label);
ylabel (ylab);
title(Title,'Interpreter', 'none')
%
hold on;
ngroups = size(Mean, 2);
%nbars = size(Mean, 2);
% Calculating the width for each bar group

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange

% Calculate center of each bar
x = (1:ngroups) ;%- groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
errorbar(x, Mean(:),STD(:), 'k', 'linestyle', 'none');

end

%
%     Trials=(RESULTSstruct.(py{d}).totalall.trial)';
%     TrialsCat=categorical(RESULTSstruct.(py{d}).totalall.trial)';
%
%     Implant=(RESULTSstruct.(py{d}).totalall.implant)';
%     ImplantCat=categorical(RESULTSstruct.(py{d}).totalall.implant)';
%
%     StepCat=categorical(RESULTSstruct.(py{d}).totalall.step)';
%     LoadCat=categorical(RESULTSstruct.(py{d}).totalall.load)';
%
%     D=RESULTSstruct.(py{d}).totalall.data';
%     N=RESULTSstruct.(py{d}).totalall.normdata';
%     N2=RESULTSstruct.(py{d}).totalall.normdata2';
%
%     T = table(D,N,Trials,Implant,StepCat,LoadCat);
%     Tstruct = table(D,N,TrialsCat,ImplantCat,StepCat,LoadCat);
%
%     AP=(Tstruct.LoadCat=='A' | Tstruct.LoadCat=='P');
%     TstructAP=Tstruct(AP,:);
%
%     IE=(Tstruct.LoadCat=='I' | Tstruct.LoadCat=='E');
%     TstructIE=Tstruct(IE,:);
%
%     I=(Tstruct.LoadCat=='I') ;
%     TstructI=Tstruct(I,:);
%
%     E=(Tstruct.LoadCat=='E') ;
%     TstructE=Tstruct(E,:);




function figuremaker( group1,  group1label,  Title, ylab, colors)

ax = gca;
V=axis;
V(4)=V(4)*1.05;
axis(V) 
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,1:length(group1));
% Naming each of the bar groups
xticklabels(ax,group1);
xlabel (group1label,'FontSize', 12);
ylabel (ylab,'FontSize', 12);
title(Title,'Interpreter', 'none','FontSize', 14)
% Creating a legend and placing it outside the bar plot
c=colors;

L1=plot(NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,'LineWidth',2);
for i=1:6
    L1(i).Color=c(i,:);
end


[leg] = legend(L1,  {'0.5','0.4', '0.3','0.2','0.1','0.6 (Patient)'},'FontSize', 11,'Location','bestoutside');
title(leg,'Conformity ratio','FontSize', 11)
leg.Title.Visible = 'on';
% lg = legend(group2,'AutoUpdate','off');
% lg.Location = 'BestOutside';
% lg.Orientation = 'Horizontal';
end

