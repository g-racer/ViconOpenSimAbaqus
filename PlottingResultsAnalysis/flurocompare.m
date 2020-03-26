load('fluro_results.mat','KNEEKIN')
types={'lungef1', 'stairupf1', 'stepupf2'};

Path = ['G:/My Drive/Thesis/ACL_Pilot_data/Processing/'];
fluropath=[Path 'jw/FLURO/'];
load([fluropath 'steptimes.mat']);


headers={'tf_FE'    'tf_VV'   'tf_IE'   'tf_ML'   'tf_AP'   'tf_SI' };

for k=1:length(types)
    trial=types{k};
    t=[];
    tf_FE=[];
    tf_VV=[];
    tf_IE=[];
    tf_ML=[];
    tf_AP=[];
    tf_SI=[];
c=1; %counter
    steptimes=fieldnames(time.(trial));
    steps=fieldnames(KNEEKIN.(trial));
    switch trial
        case 'lungef1'; no=5;
        case 'stairupf1'; no=1;
        case 'stepupf2'; no=4;
    end
    
    
    for p=1:no
        if strncmp(trial,'stair',5)==1
            csvfile=['jw_' trial '_fluoro_kinematics.csv'];
            
        else
            
            csvfile=['jw_' trial '_fluoro_kinematics_segment' num2str(p) '.csv'];
        end
        kin=csvread([fluropath csvfile],1,0);
        t=[t; kin(:,1)];

        tf_AP=[tf_AP; kin(:,3)];
        tf_SI=[tf_SI; kin(:,4)];
        tf_ML=[tf_ML; kin(:,5)];
        tf_VV=[tf_VV; kin(:,6)];
        tf_IE=[tf_IE; kin(:,7)];
        tf_FE=[tf_FE; kin(:,8)];
    end
    
    for s=1:length(steptimes)
        T=time.(trial).(steptimes{s});
        [V, idx] = min(abs(t - T));
        if min(abs(t - T))>=0.1
            disp([trial steptimes{s} ' time missing'])
        else
            
            KIN.Fluro.(trial).(steptimes{s}).t(1,1)=t(idx);
            %change in position from initial
            KIN.Fluro.(trial).(steptimes{s}).tf_AP(1,1)=tf_AP(idx)-tf_AP(1);
            KIN.Fluro.(trial).(steptimes{s}).tf_SI(1,1)=tf_SI(idx)-tf_SI(1);
            KIN.Fluro.(trial).(steptimes{s}).tf_ML(1,1)=tf_ML(idx)-tf_ML(1);
            KIN.Fluro.(trial).(steptimes{s}).tf_VV(1,1)=-tf_VV(idx);%-tf_VV(1);
            KIN.Fluro.(trial).(steptimes{s}).tf_IE(1,1)= -tf_IE(idx);%-tf_IE(1);
            KIN.Fluro.(trial).(steptimes{s}).tf_FE=-tf_FE(idx);%-tf_FE(1);
      
            %for i=[6 ] %baseline, patient
%                 try
i=6;
count=1;
                for p=1:3 %convert to degrees
                   try
                    KIN.(['implant' num2str(i)]).(trial).(steptimes{s}).(headers{p})(1,1+i)=(180/pi)*KNEEKIN.(trial).(steps{s}).(headers{p}).(['implant' num2str(i)]);
                    y=KIN.(['implant' num2str(i)]).(trial).(steptimes{s}).(headers{p})(1,1+i)...
                        -KIN.Fluro.(trial).(steptimes{s}).(headers{p})(1,1);
                    KIN.(['implant' num2str(i)]).RMSerror.(headers{p})(1,c)= y.^2; %difference squared
                   dif.(['implant' num2str(i)]).(trial)(s,p)=y;
                    count=count+1;
                   catch
%                        i=7;
%                           KIN.(['implant' num2str(i)]).(trial).(steptimes{s}).(headers{p})(1,1+i)=(180/pi)*KNEEKIN.(trial).(steps{s}).(headers{p}).(['implant' num2str(i)]);
%                     y=KIN.(['implant' num2str(i)]).(trial).(steptimes{s}).(headers{p})(1,1+i)-KIN.Fluro.(trial).(steptimes{s}).(headers{p})(1,1);
%                     KIN.(['implant' num2str(i)]).RMSerror.(headers{p})(1,c)= y.^2; %difference squared
%                     STD.(trial)(p,1)=y;
                   end
                end
                for p=4:6
                    KIN.(['implant' num2str(i)]).(trial).(steptimes{s}).(headers{p})(1,1+i)=KNEEKIN.(trial).(steps{s}).(headers{p}).(['implant' num2str(i)]);
                    y= KIN.(['implant' num2str(i)]).(trial).(steptimes{s}).(headers{p})(1,1+i)-KIN.Fluro.(trial).(steptimes{s}).(headers{p})(1,1);
                    KIN.(['implant' num2str(i)]).RMSerror.(headers{p})(1,c)=y.^2;
                  dif.(['implant' num2str(i)]).(trial)(s,p)=y;
                    count=count+1;
                end
                c=c+1;
%                 catch
%                 end
           % end
        end
        
        
    end
end

for i=[0 6] 
for p=1:length(headers)
KIN.(['implant' num2str(i)]).RMS.(headers{p})=sqrt(mean(KIN.(['implant' num2str(i)]).RMSerror.(headers{p})(1,:)));
KIN.(['implant' num2str(i)]).RMS.(headers{p})=sqrt(mean(KIN.(['implant' num2str(i)]).RMSerror.(headers{p})(1,:)));

rms.(['implant' num2str(i)])=struct2table(KIN.(['implant' num2str(i)]).RMS);

end
end







