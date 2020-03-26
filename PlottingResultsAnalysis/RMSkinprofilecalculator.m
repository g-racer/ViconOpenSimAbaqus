%RMS kin profile
clear all
close all

RMSkin=table;
types={'lunge', 'gait', 'staird', 'stairup', 'stepup'};
impnos=[1 2 3 4 5 6 ];
counter=1;
original_dir=pwd;
names=[];

for o=1:length(impnos)
    num=num2str(impnos(o));
    
    path=['jw' filesep 'RESULTS_IMPLANT' num];
    files = dir([path ]);
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    n=size(subFolders,1);
    impno=str2double({subFolders(:).name});
    impno(isnan(impno)) = [];
    impno=sort(impno);
    folderno=num2str(impno(end));
    
    imppath=[path filesep folderno];
    files = dir([imppath ]);
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    
    fnames={subFolders(:).name};
    fnames(contains(fnames(:, :), 'initial')) = [];
    fnames(contains(fnames(:, :), '.')) = [];
    fnames(contains(fnames(:, :), 'lungef1step2')) = [];
    
    apcount=1;
    iecount=1;
    for k=1:length(fnames)
        
        name=fnames{k};
        if contains(name, 'stepupf2step2')==1
            type=name(end);
            if strcmp(type, 'A')==1 || strcmp(type, 'P')==1
                file=[imppath filesep name filesep 'diag-ap.txt'];
                
                Data=dlmread(file,'');
                
                if        abs(Data(end,4))>3
                    names=[names, name];
                    continue
                else
                    %                     errorAP(apcount)=Data(end,4);
                    
                    error.AP.(['implant' num]).(name)=Data(end,4);
                    apcount=apcount+1;
                    %                     catch
                    %                         disp([name ' is empty'])
                    %                 end
                end
            else
                
                file=[imppath filesep name filesep 'diag-ie.txt'];
                Data=dlmread(file,'');
                
                if        abs(Data(end,4))>5
                    names=[names, name];
                    continue
                else
                    %                 errorIE(iecount)=Data(end,4);
                    
                    iecount=iecount+1;
                    error.IE.(['implant' num]).(name)=Data(end,4);
                    
                    %             catch
                    %                 disp([name ' is empty'])
                    %
                    %             end
                    
                    %                         %Data=D.data;
                    %                         if isempty(Data)==1
                    %                             disp([name ' is empty'])
                    %                             continue
                    %                         else
                    
                end
            end
            
            stepupf2step2.(name(end))(:,1)=Data(1:10:end,1);
            stepupf2step2.(name(end))(:,2)=Data(1:10:end,2);
            stepupf2step2.(name(end))(:,impnos(o)+2)=Data(1:10:end,3);
        end
        
    end
end



for o=1:length(impnos)
    num=num2str(impnos(o));
    try
        error.AP.(['implant' num])=rmfield(error.AP.(['implant' num]),'ngait_og1step3_A');
        error.IE.(['implant' num])=rmfield(error.IE.(['implant' num]),'ngait_og1step3_E');
        
    catch
    end
    targetRMS.AP.(['implant' num])=rms(struct2array(error.AP.(['implant' num])));
    
    targetRMS.IE.(['implant' num])=rms(struct2array(error.IE.(['implant' num])));
    
end

targetRMS.AP.all=rms(struct2array(targetRMS.AP));
targetRMS.IE.all=rms(struct2array(targetRMS.IE));

%%plots
         Colors=     [.34 .63 .83;
    .6 .4 .8;
    .3 .9 .3;
    0.88, 0.52, 0.35;
    1 .94 .0;
      .18 .52 .49]
         
%          [  0, 0.4470, 0.7410;
%         0.9290, 0.6940, 0.1250;
%         0.8500, 0.3250, 0.0980;
%         0.4660, 0.6740, 0.1880;
%         0.4940, 0.1840, 0.5560  ];


    
dir=['A' 'P'];
    figure
    hold on
    for i=1:2

    time=stepupf2step2.(dir(i))(:,1);

    plot(time, stepupf2step2.(dir(i))(:,2),'r', 'linewidth',2); %taregt
     for k=1:6
        plot(time, stepupf2step2.(dir(i))(:,k+2),'color', [Colors(k,:)], 'linewidth',1,'linestyle','--'); %taregt
     end
     
    end

    L1=plot(NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,'LineWidth',2);
L1(1).Color=[1 0 0];
    for s=2:7
    L1(s).Color=Colors(s-1,:);
       L1(s).LineStyle='--';
end
[leg] = legend(L1,  {'Target' 'CR=0.5','CR=0.4', 'CR=0.3','CR=0.2','CR=0.1','CR=0.6(Patient)'},'FontSize', 11,'Location','northwest');
% title(leg,'Legend','FontSize', 11)
% leg.Title.Visible = 'off';
xlabel ('Time (s)','FontSize', 12);
ylabel ('AP Displacement (mm)','FontSize', 12);

    
dir=['I' 'E'];
    figure
    hold on
    for i=1:2

    time=stepupf2step2.(dir(i))(:,1);

    plot(time, stepupf2step2.(dir(i))(:,2),'r', 'linewidth',2); %taregt
     for k=1:6
        plot(time, stepupf2step2.(dir(i))(:,k+2),'color', [Colors(k,:)], 'linewidth',1,'linestyle','--'); %taregt
     end
     ax=gca;
     ax.YAxisLocation='right'
     
    end
%       L1=plot(NaN,NaN,'LineWidth',2);
% L1(1).Color=[1 0 0];
%     [leg] = legend(L1,  {'Target'},'FontSize', 11,'Location','northwest');
% leg.Title.Visible = 'on';
% xlabel ('Time (s)','FontSize', 12);
ylabel ('IE rotation (deg)','FontSize', 12);

