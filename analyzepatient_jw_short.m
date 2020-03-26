
%%

patient = 'jw';


trial ='dummyboth';
name = [patient trial];

% %trial ='stairdownf2';
% trial ='stairupf1';
% trial ='stepupf2';
% trial ='lungef1';
% initial_time=0;
% final_time=1;
%path='G:/My Drive/Thesis/ACL_Pilot_data/Processing/jw/INPUTS';
%path = ['C:/Users/GRACEMCCONNOCHIE/Google_Drive/Thesis/ACL_Pilot_data/Processing/'];
Path = ['G:/My Drive/Thesis/ACL_Pilot_data/Processing/'];

cd(Path)
addpath(genpath(Path))
trials ={ 'ngait_og1', 'stairdownf2', 'stairupf1', 'stepupf2', 'lungef1'}; %

import org.opensim.modeling.* % % Import OpenSim Libraries
version = 'V6'; %Changes name of HM GS Kin Axis File
threshold_m = 0; %Minimum threshold of muscle forces
threshmax = 1500; %Maximum threshold of muscle forces
musc_multi = 1; %Multiplier of muscle forces
MPFL = '-8';
LPFL = '-5';
vs = 7843; %Call number for HM Tibia superior node
vinf = 8225; %Call number for HM Tibia inferior node

trans_path = [Path patient '/T-Matrix/']; % fill with path
scale_path = [Path patient '/Scaling/'];
trc_path= [Path patient '/MarkerData/'];
ik_path =[Path patient '/IKResults/'];
state_path =[Path patient '/StatesResults/'];
musc_amp_path= [Path patient '/LtModel/force_amp_profiles/'];
force_path= [Path patient '/ForceData/'];
ct_path = [Path patient '/LtModel/'];
STO_path=[Path patient '/STOResults/' ];
totalforce=struct;
%%
% fid = fopen([path patient '/runtimes.txt'], 'w+');
%Model
% for k=1:length(trials)
%     trial=trials{k};
name = [patient trials{k}]; %used as name in analysis transform, make sure that works
disp(['Setup for' name])
if strncmp(trial,'stair',5)==1 || strncmp(trial,'squat',5)==1 ...
        || strncmp(trial,'step',4)==1 || strncmp(trial,'chair',5)==1 || strncmp(trial,'lunge',5)==1
    
    model_file = [scale_path patient '_Model_fluro.osim'];
elseif  strncmp(trial,'dum',3)==1
    model_file = [scale_path patient '_Model_dummy.osim'];
else
    model_file = [scale_path patient '_Model.osim'];
end
model = Model([model_file]);
frame = 1;
chframe = num2str(frame);
ct_path = [Path patient '/LtModel/'];

addpath(genpath([Path '/MATLAB functions']))

%%
disp('test0');
%creates necessary setup files

force_file = [force_path trial '.mot'];
kin_F_file = [ik_path trial '_F_ik.mot'];
STOfile=[STO_path trial '_StaticOptimization_force.sto'];
state_file = [state_path trial '_StatesReporter_states.sto'];
[ initial_time, final_time ] = gettrialtime(kin_F_file,force_file,patient, trial);
[time ] = getstoptimes(Path, patient,trial, initial_time, final_time);

% fprintf(fid,['Trial: ' name '\n']);
% MyFieldNames = fieldnames(time);
% for i=1:4
%     p= getfield(time,MyFieldNames{i});
% fprintf(fid,'%1.4f\n,', p);
% end
% fclose(fid);

% fill with output file of IK (.mot)
if exist(kin_F_file,'file') ~= 2
    disp('Running IK full cycle analysis')
    setupAndRunIKSingleCycle( patient,Path, trial)
end
% time  = getstoptimes(path, patient,trial, initial_time, final_time)

%run cropped inverse kinematics
kin_file = [ik_path trial '_ik.mot']; % fill with output file of IK (.mot)
if exist(kin_file,'file') ~= 2
    disp('Running IK cropped cycle analysis')
    setupAndRunIKSingleCycle( patient,Path, trial,initial_time, final_time )
end


% output of the StateReporter Analyze tool
if exist(state_file,'file') ~= 2
    disp('Running States analysis')
    setupAndRunStatesSingleCycle( patient,Path, trial, initial_time, final_time)
end


if exist(STOfile,'file') ~= 2
    disp('Running STO analysis')
    setupAndRunSTOSingleCycle(patient, Path, trial, initial_time, final_time)
end





%%
disp('test1')
%Transformation Matrices
%Creates State Reporter File and Calculates the Tranformation Matrices

%Requires Model's osim, trc, mot Files, and generic SetUp SR file
%Input patient and surgery
%Ensure 'Body Name' is correct for all T
%Run in Double Knee folder


%Determine time range for one gait cycle
%[ initial_time, final_time ] = gettrialtime(kin_file,force_file, trial);

% Get trc data to determine time range
% trc_file = [trc_path trial '.trc'];
% markerData = MarkerData(trc_file);
% % Get initial and intial time
% initial_time = markerData.getStartFrameTime();
% final_time = markerData.getLastFrameTime();



if exist([trans_path name '.mat'],'file')==2
    load([trans_path name '.mat']);
else
    
    %Calculate the transformation matricies for each body and each frame
    %using States file
    %bodies are different names from bones
    %[T_Matrix] = analysis_transform(model, storage, fcn, name, sr_path)
    %fcn= Transform_body(model,state,body_name)
    %analysis_transform(model, storage, fcn, name, state_path)
    
    [T_pel] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'pelvis'),trial,state_path);
    [T_fem] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'femur'),trial,state_path);
    [T_tib] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'tibia'),trial,state_path);
    [T_pat] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'patella'),trial,state_path);
    [T_tal] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'talus'),trial,state_path);
    [T_cal] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'calcaneus'),trial,state_path);
    [T_toe] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'toes'),trial,state_path);
    [T_FemComp] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'femoral_component'),trial,state_path);
    [T_TTray] = analysis_transform(model_file,state_file,...
        @(model, state) Transform_body(model,state,'tibial_tray'),trial,state_path);
    TM=struct;
    TM.T_pel = T_pel.pelvis(:,:,:);
    TM.T_fem = T_fem.('femur')(:,:,:);
    TM.T_tib = T_tib.('tibia')(:,:,:);
    TM.T_pat = T_pat.('patella')(:,:,:);
    TM.T_tal = T_tal.('talus')(:,:,:);
    TM.T_cal = T_cal.('calcaneus')(:,:,:);
    TM.T_toe = T_toe.('toes')(:,:,:);
    TM.T_FemComp = T_FemComp.('femoral_component')(:,:,:);
    
    TM.T_TTray = T_TTray.('tibial_tray')(:,:,:);
    % TM.T_tib(:,[1],:)=TM.T_TTray(:,[1],:)
    %TM.T_tib(:,[2 4],:)=TM.T_TTray(:,[2 4],:)*-1;
    
    save([trans_path name '.mat'],'TM'); %'-struct',
  
    %load([trans_path name '.mat'],'T*');
end




%%
disp('test2')
%Bone Coordinate Transformation
%Creates Bone inp Files
%Requires Generic Bone inp Files
%Input frame
%Have Run TranformationMatrices.m
%Ensure if elseif statement for bone and body is correct


%creates a new folder for the frame of interest to put input files states in this frame
if exist(([ct_path name  'Frame' chframe]),'dir')==0;
    mkdir([ct_path name  'Frame' chframe]);
end
%load transformed node coordinates
if exist([trans_path name 'frame' num2str(frame) 'coord.mat'],'file')==2
    load([trans_path name 'frame' num2str(frame) 'coord.mat'])
    load([trans_path name 'frame' num2str(frame) 'GScoord.mat']);
    load([trans_path name 'frame' num2str(frame) 'Refcoord.mat']);
    
else
    
    Nodescord=struct;
    GSnodes=struct;
    REFnodes=struct;
    for h = [3 11  ]
        if h == 1
            body = 'pelvis';
            bone = 'pelvis';
        elseif h == 2
            body = 'femur';
            bone = 'femur';
        elseif h == 3
            body = 'patella';
            bone = 'patella';
        elseif h == 4
            body = 'tibia';
            bone = 'tibia';
        elseif h == 6
            body = 'tibia';
            bone = 'fibula';
        elseif h == 5
            body = 'talus';
            bone = 'talus';
        elseif h == 7
            body = 'calcn';
            bone = 'foot';
        elseif h == 8
            body = 'toes';
            bone = 'bofoot';
        elseif h == 9
            body = 'femoral_component';
            bone = 'femoral_component';
        elseif h == 10
            body = 'tibial_tray';
            bone = 'tibial_tray';
        elseif h == 11
            body = 'patella';
            bone = 'dome';
        else
            break
        end
        
        switch body %select correct transformation matix
            case 'pelvis'; T = TM.T_pel(:,:,frame);
            case 'femur'; T = TM.T_fem(:,:,frame);
            case 'tibia'; T = TM.T_tib(:,:,frame);
            case 'patella'; T = TM.T_pat(:,:,frame);
            case 'dome'; T = TM.T_pat(:,:,frame);
            case 'talus'; T = TM.T_tal(:,:,frame);
            case 'calcn'; T = TM.T_cal(:,:,frame);
            case 'toes'; T = TM.T_toe(:,:,frame);
            case 'femoral_component'; T = TM.T_FemComp(:,:,frame);
            case 'tibial_tray'; T = TM.T_TTray(:,:,frame);
        end
        
        %Uses boneread function to write input files and transform node
        %coordinates.
        if strncmp(bone, 'femur',5)==1 || strncmp(bone, 'patella',7)==1 || ...
                strncmp(bone, 'talus',5)==1 || strncmp(bone, 'pelvis',6)==1 ||...
                strcmp(bone, 'tibia')==1
            [nodes, gsnodes, refnodes]=bonereadjw(patient,trial,bone,Path,T,frame);
            Nodescord.(bone)=nodes;
            GSnodes.(bone)= gsnodes;
            REFnodes.(bone)= refnodes;
            %             if strncmp(bone,'fem',3)==1 || strncmp(bone,'pelvis',5)==1
            %                 GSnodes.(bone)(:,3)= gsnodes(:,3)-12;
            %                 REFnodes.(bone)= refnodes;
            %             elseif strncmp(bone,'patella',7)==1
            %                 GSnodes.(bone)(:,3)= gsnodes(:,3)-4;
            %                 GSnodes.(bone)(:,2)= gsnodes(:,2)+3;
            %                 REFnodes.(bone)= refnodes;
            %             end
        else
            nodes=bonereadjw(patient,trial,bone,Path,T,frame);
            Nodescord.(bone)=nodes;
        end
    end
    % T = TM.T_fem(:,:,frame);
    % ligreadjw(patient,trial,path,T,frame);
    % Nodescord.('lig')=nodes;
    
    save([trans_path name 'frame' num2str(frame) 'coord.mat'],'Nodescord');
    save([trans_path name 'frame' num2str(frame) 'GScoord.mat'],'GSnodes');
    save([trans_path name 'frame' num2str(frame) 'Refcoord.mat'],'REFnodes');
    
end
% %%
% disp('test3 -muscleframecoordinates')
% %create muscle connector elements for frame/s
% frame=1;
% muscleframecoordinates(path,patient,trial,frame);



%%
disp('test4 -calc kinematics')
%Transform coordinates to global coordinate system of each frame; combine into one matrix
%to create kinematic amplitude card for gait cycle
gs = zeros(20,4,size(TM.T_tib,3));
Bones ={    'femur' 'tibia' 'patella' 'pelvis',  'talus'  }; %

for i=1:5
    bone=Bones{i};
    switch i %Select Bone Coordinates and Transformation Matrix
        case 1; T = TM.T_fem;
        case 2; T = TM.T_tib;
        case 3;   T = TM.T_pat;
        case 4;  T = TM.T_pel;
        case 5;   T=TM.T_tal;
    end
    
    for k = 1:size(T,3) %num frames
        for j=1:4
            %for each GS node, and for each frame, multiply coordinates by transformation matrix
            GSTrans.(bone)(j,:,k)=T(:,:,k)*([1 GSnodes.(bone)(j,2:4)])';
            gs(4*i+j-4,:,k)= GSTrans.(bone)(j,:,k);%T(:,:,k)*([1 GSnodes.(bone)(j,2:4)])';
        end
    end
end

gs_Gen = zeros(4,4,size(T,3)); gs_Gen(:,1,:) = 1; %Blank GS Matrix
gsR_Fem = gs_Gen; gsR_Fem = gs(1:4,:,:);    %Right Femur GS Matrix
gsR_Tib = gs_Gen; gsR_Tib = gs(5:8,:,:);    %Right Tibia GS Matrix
gsR_Pat = gs_Gen; gsR_Pat = gs(9:12,:,:);   %Right Patella GS Matrix
gs_Pel = gs_Gen; gs_Pel = gs(13:16,:,:);    %Pelvis GS Matrix, Right node treated as Medial
gsR_Tal = gs_Gen; gsR_Tal = gs(17:20,:,:);  %Right Talus GS Matrix


% figure
% for i = 1:size(T,3)
%     for j=1:4
%         hold on
%         plot3([gsR_Fem(1,4,i),gsR_Fem(2,4,i)],[gsR_Fem(1,2,i),gsR_Fem(2,2,i)],[gsR_Fem(1,3,i),gsR_Fem(2,3,i)],'r')
%             plot3([gsR_Fem(3,4,i),gsR_Fem(4,4,i)],[gsR_Fem(3,2,i),gsR_Fem(4,2,i)],[gsR_Fem(3,3,i),gsR_Fem(4,3,i)],'r')
%         plot3([gsR_Tib(1,4,i),gsR_Tib(2,4,i)],[gsR_Tib(1,2,i),gsR_Tib(2,2,i)],[gsR_Tib(1,3,i),gsR_Tib(2,3,i)],'g')
%             plot3([gsR_Tib(3,4,i),gsR_Tib(4,4,i)],[gsR_Tib(3,2,i),gsR_Tib(4,2,i)],[gsR_Tib(3,3,i),gsR_Tib(4,3,i)],'b')
% %         scatter3(gsR_Fem(j,4,i),gsR_Fem(j,2,i),gsR_Fem(j,3,i),'r')
% %     scatter3(gsR_Tib(j,4,i),gsR_Tib(j,2,i),gsR_Tib(j,3,i),'b' );
%     end
%     drawnow
%     pause(0.2)
% % end
% figure
% hold on
% plot3([gsR_Fem(1,4,1),gsR_Fem(2,4,1)],[gsR_Fem(1,2,1),gsR_Fem(2,2,1)],[gsR_Fem(1,3,1),gsR_Fem(2,3,1)],'r') %ML
%             plot3([gsR_Fem(3,4,1),gsR_Fem(4,4,1)],[gsR_Fem(3,2,1),gsR_Fem(4,2,1)],[gsR_Fem(3,3,1),gsR_Fem(4,3,1)],'b') %IE
%         plot3([gsR_Tib(1,4,1),gsR_Tib(2,4,1)],[gsR_Tib(1,2,1),gsR_Tib(2,2,1)],[gsR_Tib(1,3,1),gsR_Tib(2,3,1)],'g')
%             plot3([gsR_Tib(3,4,1),gsR_Tib(4,4,1)],[gsR_Tib(3,2,1),gsR_Tib(4,2,1)],[gsR_Tib(3,3,1),gsR_Tib(4,3,1)],'m')

kin = zeros(size(T,3),24);
inter_rtf1=zeros(size(T,3),3);
inter_rtf2=zeros(size(T,3),3);
inter_rpf1=zeros(size(T,3),3);
inter_rpf2=zeros(size(T,3),3);
inter_rhip1=zeros(size(T,3),3);
inter_rhip2=zeros(size(T,3),3);
inter_rank1=zeros(size(T,3),3);
inter_rank2=zeros(size(T,3),3);
for ii = 1:size(T,3)
    %      [kin(ii,1:6), inter_rtf1(ii,:), inter_rtf2(ii,:)] = GroodSuntay(gsR_Fem(:,:,ii), gsR_Tib(:,:,ii), 'Inf','Sup','top'); %changed from ML-fix femur
    %      [kin(ii,7:12), inter_rpf1(ii,:), inter_rpf2(ii,:)] = GroodSuntay(gsR_Fem(:,:,ii), gsR_Pat(:,:,ii) ,'Inf','Inf','top'); %fix femur
    %      [kin(ii,13:18), inter_rhip1(ii,:), inter_rhip2(ii,:)] = GroodSuntay( gsR_Fem(:,:,ii),gs_Pel(:,:,ii),'Sup','Med','BoneB'); %move femur wrt pelvis
    %    [kin(ii,19:24), inter_rank1(ii,:), inter_rank2(ii,:)] = GroodSuntay(gsR_Tib(:,:,ii), gsR_Tal(:,:,ii),'Inf','Sup','top');
    %original
    [kin(ii,1:6), inter_rtf1(ii,:), inter_rtf2(ii,:)] = GroodSuntay(gsR_Tib(:,:,ii),gsR_Fem(:,:,ii),'Sup','Inf','BoneB'); %,X_si_A(ii,:),X_ml_B(ii,:),X_ap_F(ii,:)
    [kin(ii,7:12), inter_rpf1(ii,:), inter_rpf2(ii,:)] = GroodSuntay( gsR_Pat(:,:,ii), gsR_Fem(:,:,ii),  'Sup','Inf','BoneB');
    [kin(ii,13:18), inter_rhip1(ii,:), inter_rhip2(ii,:)] = GroodSuntay( gsR_Fem(:,:,ii), gs_Pel(:,:,ii),'Sup','Med','BoneB'); %changes to bone B
    [kin(ii,19:24), inter_rank1(ii,:), inter_rank2(ii,:)] = GroodSuntay(  gsR_Tib(:,:,ii), gsR_Tal(:,:,ii),  'Inf', 'Sup', 'BoneA');
    % axistib(1:2,1:3,ii)=[gsR_Tib(4,2:4,ii); gsR_Tib(4,2:4,ii)-400*X_si_A(ii,:)];
    % axisfem(1:2,1:3,ii)=[gsR_Fem(3,2:4,ii); gsR_Fem(3,2:4,ii)+10*X_ml_B(ii,:)];
    % plot3([axistib(1,3,ii),axistib(2,3,ii)],[axistib(1,1,ii),axistib(2,1,ii)],[axistib(1,2,ii),axistib(2,2,ii)],'b')
    % plot3([axisfem(1,3,ii),axisfem(2,3,ii)],[axisfem(1,1,ii),axisfem(2,1,ii)],[axisfem(1,2,ii),axisfem(2,2,ii)],'r')
    % %
    % plot3([gsR_Fem(1,4,ii),gsR_Fem(2,4,ii)],[gsR_Fem(1,2,ii),gsR_Fem(2,2,ii)],[gsR_Fem(1,3,ii),gsR_Fem(2,3,ii)],'m')
    %             plot3([gsR_Fem(3,4,ii),gsR_Fem(4,4,ii)],[gsR_Fem(3,2,ii),gsR_Fem(4,2,ii)],[gsR_Fem(3,3,ii),gsR_Fem(4,3,ii)],'m')
    %         plot3([gsR_Tib(1,4,ii),gsR_Tib(2,4,ii)],[gsR_Tib(1,2,ii),gsR_Tib(2,2,ii)],[gsR_Tib(1,3,ii),gsR_Tib(2,3,ii)],'g')
    %             plot3([gsR_Tib(3,4,ii),gsR_Tib(4,4,ii)],[gsR_Tib(3,2,ii),gsR_Tib(4,2,ii)],[gsR_Tib(3,3,ii),gsR_Tib(4,3,ii)],'g')
    %              quiver3(gsR_Fem(3,4,1),gsR_Fem(3,2,1),gsR_Fem(3,3,1),inter_rtf2(ii,3),inter_rtf2(ii,1),inter_rtf2(ii,2),'r')
    % quiver3(gsR_Tib(3,4,1),gsR_Tib(3,2,1),gsR_Tib(3,3,1),inter_rtf1(ii,3),inter_rtf1(ii,1),inter_rtf1(ii,2),'b')
    %
    %
    % drawnow
    %     pause(0.1)
end
kin(:,[1:3 7:9 13:15 19:21]) = unwrap(kin(:,[1:3 7:9 13:15 19:21]),pi/2); %fix jump in angle
t = linspace(initial_time,final_time, size(kin,1));
t=t-initial_time;

%%%%%
%kin=kin-kin(1,:); %make relative??

close all
% figure
% plot(t,kin(:,1:3))
% legend('ML', 'AP' ,'SI')
% title('knee angle')
% figure
% plot(t,kin(:,4:6)-kin(1,4:6))
% legend('ML', 'AP' ,'SI')
% title('knee disp')
%
% %  for i=1:3
% %      figure
% %      plot(t,kin(:,i+2))
% %  end
% %
%
% %
% figure
% plot(t,kin(:,7:9))
% legend('ML', 'AP' ,'SI')
% title('PF angle')
%
% figure
% plot(t,kin(:,10:12)-kin(1,10:12))
% legend('ML', 'AP' ,'SI')
% title('PF disp')
%
%
% figure
% plot(t,kin(:,13:15))
% legend('ML', 'AP' ,'SI')
% title('hip angle')
%
% figure
% plot(t,kin(:,16:18))
% legend('ML', 'AP' ,'SI')
% title('hip disp')
%
%
% figure
% plot(t,kin(:,19:21))
% legend('ML', 'AP' ,'SI')
% title('ankle angle')
%
% figure
% plot(t,kin(:,22:24)-kin(1,22:24))
% legend('ML', 'AP' ,'SI')
% title('ankle disp')
%
% % close all



% quiver3(0,0,0,inter_rhip1(1,1),inter_rhip1(1,2),inter_rhip1(1,3))
% quiver3(0,0,0,inter_rhip2(1,1),inter_rhip2(1,2),inter_rhip2(1,3))

%%
disp('test5 -Write Kinematic Amplitude Card for full cycle')
%Write Kinematic Amplitude Card for full cycle
if exist([Path patient '\LtModel\amp_profiles\' trial '\']) ==0;
    mkdir([Path patient '\LtModel\amp_profiles\' trial '\']);
end
f = kin(:,[1:6 13:24]); %ignore patella %kinematics
T = linspace(initial_time,final_time, size(f,1));
t = (T-initial_time); % set to start at 0s%=0.025/2
%SCALED
%t=t*10;
% t = initial_time:time_increment:final_time;
% t = (t-initial_time)*3; % set to start at 0s%=0.025/2

% f = kin;
% t = 0:time_increment:final_time;
% t = t*3;
% f(:,13:18)=f(:,13:18)*-1;
fout1=fopen([Path patient '\LtModel\amp_profiles\' trial '\' trial '-KIN.inp'],'w'); %

amp_list = {'AMP_RTF_FE_LOAD'; 'AMP_RTF_VV_LOAD'; 'AMP_RTF_IE_LOAD'; 'AMP_RTF_ML_LOAD'; 'AMP_RTF_AP_LOAD'; 'AMP_RTF_SI_LOAD'; ...
    %  'AMP_RPF_FE_LOAD'; 'AMP_RPF_VV_LOAD'; 'AMP_RPF_IE_LOAD'; 'AMP_RPF_ML_LOAD'; 'AMP_RPF_AP_LOAD'; 'AMP_RPF_SI_LOAD'; ...
    'AMP_RHip_FE_LOAD'; 'AMP_RHip_VV_LOAD'; 'AMP_RHip_IE_LOAD'; 'AMP_RHip_ML_LOAD'; 'AMP_RHip_AP_LOAD'; 'AMP_RHip_SI_LOAD'; ...
    'AMP_RAnk_FE_LOAD'; 'AMP_RAnk_VV_LOAD'; 'AMP_RAnk_IE_LOAD'; 'AMP_RAnk_ML_LOAD'; 'AMP_RAnk_AP_LOAD'; 'AMP_RAnk_SI_LOAD'};

for var = 1:size(amp_list,1)
    
    NAME = amp_list{var};
    fprintf(fout1,['*AMPLITUDE, VALUE=RELATIVE, SMOOTH=0.1, NAME=' NAME '\n']);
    for pc=1:size(f,1) %no. frames
        if rem(pc,4)==1; fprintf(fout1,'%-1.4f,%-1.4f',[t(pc) f(pc,var)]);
        elseif rem(pc,4)==0; fprintf(fout1,',%-1.4f,%-1.4f\n',[t(pc) f(pc,var)]);
        else
            fprintf(fout1,',%-1.4f,%-1.4f',[t(pc) f(pc,var)]);
        end
    end
    if rem(pc,4)>0; fprintf(fout1,'\n**\n'); %rem =remainder
    else
        fprintf(fout1,'**\n');
    end
end

fclose(fout1);
%%
disp('test6 -Write Kinematic Amplitude Card for each step and hold step')
%Write Kinematic Amplitude Card for each step and hold step
load([ct_path 'initialstance/initialkin.mat'])
%kinematics
T = linspace(initial_time,final_time, size(kin,1));
for i=1:4 %4 steps
    num=num2str(i);
    fout1=fopen([Path patient '\LtModel\amp_profiles\' trial '\' trial '_KIN' num '.inp'],'w'); %
    t1=initial_time;
    if contains(trial,'lunge')==1
        switch i
            case 1;  t2=time.heelstrike; %kinematics
            case 2;  t2=time.peak1;
            case 3;  t2=time.maxflex;
            case 4;  t2=time.toeoff;
        end
    else
        
        switch i
            case 1;  t2=time.heelstrike; %kinematics
            case 2;  t2=time.peak1;
            case 3;  t2=time.peak2;
            case 4;  t2=time.toeoff;
        end
    end
    
    [V, idxend] = (min(abs(T - t2)));
    [V, idxstrt] =(min(abs(T - t1)));
    ind=idxstrt:idxend;
    
    f = kin(ind,[1:6 13:24]); %only want time period of interest and ignore patella
    tt=linspace(t1,t2, size(f,1));
    tt = (tt-t1);%initialize to 0 for step time
    %SCALED
    SimTime=1.75;
    SFactor = SimTime/tt(end);  % scaling factor
    tt = tt*SFactor;
    
    amp_list = {'AMP_RTF_FE_LOAD'; 'AMP_RTF_VV_LOAD'; 'AMP_RTF_IE_LOAD'; 'AMP_RTF_ML_LOAD'; 'AMP_RTF_AP_LOAD'; 'AMP_RTF_SI_LOAD'; ...
        %  'AMP_RPF_FE_LOAD'; 'AMP_RPF_VV_LOAD'; 'AMP_RPF_IE_LOAD'; 'AMP_RPF_ML_LOAD'; 'AMP_RPF_AP_LOAD'; 'AMP_RPF_SI_LOAD'; ...
        'AMP_RHip_FE_LOAD'; 'AMP_RHip_VV_LOAD'; 'AMP_RHip_IE_LOAD'; 'AMP_RHip_ML_LOAD'; 'AMP_RHip_AP_LOAD'; 'AMP_RHip_SI_LOAD'; ...
        'AMP_RAnk_FE_LOAD'; 'AMP_RAnk_VV_LOAD'; 'AMP_RAnk_IE_LOAD'; 'AMP_RAnk_ML_LOAD'; 'AMP_RAnk_AP_LOAD'; 'AMP_RAnk_SI_LOAD'};
    %print amplitude card for step
    for var = 1:size(amp_list,1)
        NAME = [amp_list{var}];
       % fprintf(fout1,['*AMPLITUDE, VALUE=RELATIVE, SMOOTH=0.2, NAME=' NAME '\n']);
                fprintf(fout1,['*AMPLITUDE, VALUE=RELATIVE, DEFINITION=SMOOTH STEP,  NAME=' NAME '\n']);

        %  initialamp.(NAME)=f(1,var);
        %         for pc=1:size(f,1) %no. frames
        %             if rem(pc,4)==1; fprintf(fout1,'%-1.4f,%-1.4f',[tt(pc) f(pc,var)]);
        %             elseif rem(pc,4)==0; fprintf(fout1,',%-1.4f,%-1.4f\n',[tt(pc) f(pc,var)]);
        %             else
        %                 fprintf(fout1,',%-1.4f,%-1.4f',[tt(pc) f(pc,var)]);
        %             end
        %         end
        %         if rem(pc,4)>0; fprintf(fout1,'\n**\n'); %rem =remainder
        %         else
        %             fprintf(fout1,'**\n');
        %         end
        fint=initialamp.(NAME);
        %   fprintf(fout1,'%1.4f, %1.4f\n',[tt(1) f(1,var)]);
        %    fprintf(fout1,'%1.4f, %1.4f\n',[0.1 f(1,var) ]);
        fprintf(fout1,'%1.4f, %1.4f\n',[0 fint]);
        % fprintf(fout1,'%1.4f, %1.4f\n',[0.1 fint ]);
        %
        %     if contains(NAME,'AMP_RTF_AP_LOAD')==1
        %         fprintf(fout1,'%1.4f, %1.4f\n',[ 1.25 0 ]);
        %          fprintf(fout1,'%1.4f, %1.4f\n',[ tt(end) 0]);
        %     else
        fprintf(fout1,'%1.4f, %1.4f\n',[ 1.7 f(end,var) ]);
        fprintf(fout1,'%1.4f, %1.4f\n',[ 2.1 f(end,var)]);
        %     end
        fprintf(fout1,'**\n');
        if var==1
            Kneeflex.(trial).(['step' num2str(i)])=f(end,var);
        end
    end
    % save([ct_path 'initialstance/initialkin.mat'],'initialamp');
    fclose(fout1);
end
save([musc_amp_path    '/Kneeflex.mat'],'Kneeflex');

%     %print kin amplitude card for AP and SI for hold
%     %%%AP LOAD
%     fout2=fopen([path patient '\LtModel\amp_profiles\' trial '\' trial '_holdAP-KIN.inp'],'w'); %
% %     NAME = ['AMP_RTF_AP_LOAD' num '_hold'];
% %     %     fprintf(fout2,['*AMPLITUDE,  DEFINITION=PERIODIC, NAME=' NAME '\n']);
% %     %     fprintf(fout2,'%-1.0f, %-1.4f, %-1.4f, %-1.4f\n',[1 2*pi 0.0 0.0 ]); %N f to Ao
% %     %     fprintf(fout2,'%-1.4f, %-1.4f\n',[0 1 ]); %A and B cos and sin?
% %
% %     fprintf(fout2,['*AMPLITUDE,  DEFINITION=SMOOTH STEP, NAME=' NAME '\n']);
% %     fprintf(fout2,'%-1.0f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f\n' ...
% %         ,[0 0 0.25 1 0.5 0 0.75 -1 ]);
% %     fprintf(fout2, '%-1.4f, %-1.4f\n', [1 0]);
% %     fprintf(fout2,'**\n');
% %
% %     fprintf(fout2,['** Name:' NAME '_load  Type: Connector force\n']);
% %     fprintf(fout2,['*CONNECTOR LOAD, OP=NEW, amplitude=' NAME '\n']);
% %     fprintf(fout2,'AXIS_R_TF_EXT_AP, ');
% %     fprintf(fout2,' %-1.0f, \t %-1.4f\n',[1 -200]); %200N AP along axis direction
%
%
%     axes_list =  ...
%         { 'R_Hip_ML', 'R_Hip_AP',  'R_Hip_IE','R_Hip_ML', 'R_Hip_AP',  'R_Hip_IE',   ...
%         'R_Ank_ML', 'R_Ank_AP','R_Ank_IE',   'R_Ank_ML', 'R_Ank_AP','R_Ank_IE',...
%         'R_TF_EXT_ML', 'R_TF_EXT_IE' };
%     for var=1:size(axes_list,2)
%
%
%         fprintf(fout2,['*CONNECTOR MOTION, OP=NEW\n' ]); %Change from load FIXED, OP=NEW, amplitude=hold
%         fprintf(fout2,['axis_' axes_list{var} ','] );
%         id=[1 2 3 7 8 9 13 ]; %angles
%         if any(var==id)==1
%             fprintf(fout2, '%-1.0f, %-1.4f\n', [4 0]  ); %rotation first f(end,var)] , \t %-1.4f
%         else
%             fprintf(fout2, '%-1.0f, %-1.4f\n', [1 0] ); %translation  f(end,var)]
%
%         end
%     end
%     fclose(fout2);  %
%
%
%     %%%%%%%%%%%%%IE
%
%     fout2=fopen([path patient '\LtModel\amp_profiles\' trial '\' trial  '_holdIE-KIN.inp'],'w'); %
% %     NAME = ['AMP_RTF_IE_LOAD' num '_hold'];
% %     %     fprintf(fout2,['*AMPLITUDE,  DEFINITION=PERIODIC, NAME=' NAME '\n']);
% %     %     fprintf(fout2,'%-1.0f, %-1.4f, %-1.4f, %-1.4f\n',[1 2*pi 0.0 0.0 ]); %N f to Ao
% %     %     fprintf(fout2,'%-1.4f, %-1.4f\n',[0 1 ]); %A and B cos and sin?
% %     fprintf(fout2,['*AMPLITUDE,  DEFINITION=SMOOTH STEP, NAME=' NAME '\n']);
% %     fprintf(fout2,'%-1.0f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f\n' ...
% %         ,[0 0 0.25 1 0.5 0 0.75 -1 ]);
% %     fprintf(fout2, '%-1.4f, %-1.4f\n', [1 0]);
% %     fprintf(fout2,'**\n');
% %
% %     fprintf(fout2,['** Name:' NAME '_load  Type: Connector force\n']);
% %     fprintf(fout2,['*CONNECTOR LOAD, OP=NEW, amplitude=' NAME '\n']);
% %     fprintf(fout2,'AXIS_R_TF_EXT_IE, ');
% %     fprintf(fout2,' %-1.0f, \t %-1.4f\n',[1 5000]); %5000Nmm IE along axis direction
%
%     axes_list =         { 'R_Hip_ML', 'R_Hip_AP',  'R_Hip_IE','R_Hip_ML', 'R_Hip_AP',  'R_Hip_IE',   ...
%         'R_Ank_ML', 'R_Ank_AP','R_Ank_IE',   'R_Ank_ML', 'R_Ank_AP','R_Ank_IE',...
%         'R_TF_EXT_ML', 'R_TF_EXT_AP' };
%     for var=1:size(axes_list,2)
%
%
%         fprintf(fout2,['*CONNECTOR MOTION, OP=NEW\n' ]); %Change from load FIXED, OP=NEW, amplitude=hold
%         fprintf(fout2,['axis_' axes_list{var} ','] );
%         id=[1 2 3 7 8 9 13 ]; %angles
%         if any(var==id)==1
%             fprintf(fout2, '%-1.0f\n', [4 ]  ); %rotation first f(end,var)] , \t %-1.4f
%         else
%             fprintf(fout2, '%-1.0f\n', [1 ] ); %translation  f(end,var)]
%
%         end
%     end
%     for var=1:size(axes_list,2)
%         fprintf(fout2,['*CONNECTOR MOTION,  OP=NEW\n' ]); %Change from load FIXED,
%         fprintf(fout2,['axis_' axes_list{var} ','] );
%         id=[1 2 3 7 8 9 13 ]; %angles
%         if any(var==id)==1
%             fprintf(fout2, '%-1.0f, %-1.4f\n', [4 0]  ); %rotation first f(end,var)] , \t %-1.4f
%         else
%             fprintf(fout2, '%-1.0f, %-1.4f\n', [1 0] ); %translation  f(end,var)]
%
%         end
%     end



%
%                  amp_list = {'AMP_RTF_IE_LOAD'; 'AMP_RTF_AP_LOAD'};
%                  for var = 1:size(amp_list,1)
%                      NAME = [amp_list{var} num '_hold'];
%                      %         if contains(amp_list{var},'RTF_IE')==1 || contains(amp_list{var},'RTF_AP')==1 %apply AP IE loads
%
%
%
%
%         if contains(amp_list{var},'RTF_AP')==1
%
%             fprintf(fout2,['** Name:' NAME '_load  Type: Connector force\n']);
% %             fprintf(fout2,['*CLoad, OP=NEW, amplitude=' NAME '\n']);
% %             fprintf(fout2,'200, ');
%                         fprintf(fout2,['*Connector Load, OP=NEW, amplitude=' NAME '\n']);
%                         fprintf(fout2,'AXIS_R_TF_EXT_AP, ');
%             fprintf(fout2,' %-1.0f, \t %-1.4f\n',[1 -200]); %500N AP along axis direction
%         else
% %             fprintf(fout2,'%-1.0f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f, %-1.4f\n' ...
% %                 ,[0 0 1 0 1.25 1 1.5 0 1.75 -1 2 0]);
%             fprintf(fout2,'**\n');
%             fprintf(fout2,['** Name:' NAME '_load  Type: Connector moment\n']);
% %             fprintf(fout2,['*CLoad, OP=NEW, amplitude=' NAME '\n']);
% %             fprintf(fout2,'200, '); %RBRN
%                         fprintf(fout2,['*Connector Load, OP=NEW, amplitude=' NAME '\n']);
%                         fprintf(fout2,'AXIS_R_TF_EXT_IE, ');
%             fprintf(fout2,' %-1.0f, \t %-1.4f\n',[4 5000]); %axis direction of rotation 5000Nmm IE
%         end
%         fprintf(fout2,'**\n');
%     end

%                 fprintf(fout2,'%-1.4f,%-1.4f, %-1.4f,%-1.4f %-1.4f,%-1.4f, %-1.4f,%-1.4f,\n',[0 0 0.25 10 0.5 0 0.75 -10 ]);
%              fprintf(fout2,' %-1.4f,%-1.4f\n',[1 0]);
%             else
%             fprintf(fout2,'%-1.4f,%-1.4f, %-1.4f,%-1.4f %-1.4f,%-1.4f, %-1.4f,%-1.4f,\n',[0 0 0.25 500 0.5 0 0.75 -500 ]);
%              fprintf(fout2,' %-1.4f,%-1.4f\n',[1 0]);
%             end
%             %              fprintf(fout2,'**\n');
%         else
%             fprintf(fout2,['*AMPLITUDE, NAME=' NAME  ', DEFINITION=SMOOTH STEP\n']); %VALUE=RELATIVE, SMOOTH=0.1, DEFINITION=SMOOTH STEP,
%             fprintf(fout2,'%-1.4f,\t%-1.4f\n%-1.4f,\t%-1.4f\n',[0 f(end,var) 1 f(end,var)]);
%             fprintf(fout2,'**\n');

%             %fix axes
%     axes_list = ...%{'R_TF_EXT_ML', 'R_TF_EXT_AP', 'R_TF_EXT_IE', 'R_TF_EXT_ML', 'R_TF_EXT_AP', 'R_TF_EXT_IE'  ...
%         { 'R_Hip_ML', 'R_Hip_AP',  'R_Hip_IE','R_Hip_ML', 'R_Hip_AP',  'R_Hip_IE',   ...
%         'R_Ank_ML', 'R_Ank_AP','R_Ank_IE',   'R_Ank_ML', 'R_Ank_AP','R_Ank_IE', };
%     for var=1:size(axes_list,2)
%         fprintf(fout2,['*CONNECTOR MOTION, FIXED, OP=NEW\n' ]); %Change from load
%         fprintf(fout2,['axis_' axes_list{var} ','] );
%         id=[1 2 3 7 8 9 13 14 15]; %angles
%         if any(var==id)==1
%             fprintf(fout2, ' %-1.0f, %-1.4f\n', [4 f(end,var+6)]  ); %rotation first f(end,var)] , \t %-1.4f
%         else
%             fprintf(fout2, ' %-1.0f, %-1.4f\n', [1 f(end,var+6)] ); %translation  f(end,var)]
%
%         end
%         fprintf(fout2,'**\n');
%     end
%     fclose(fout2);  %
% end


fclose('all');





%%
%%MUSCLE amp file
disp('test8 -MUSCLE amp file')
%
Musc_amp_list = {'patlig_lat' 'patlig_center' 'patlig_med' 'addbrev' 'addlong' 'addmagProx' 'addmagMid' 'addmagDist' 'addmagIsch'...
    'bflh' 'bfsh' 'edl' 'ehl' 'fdl' 'fhl' 'gaslat' 'gasmed' 'gem'...
    'glmax1' 'glmax2' 'glmax3' 'glmed1' 'glmed2' 'glmed3' 'glmin1' 'glmin2' 'glmin3' ...
    'grac' 'iliacus' 'pect' 'perbrev' 'perlong' 'pertert' 'piri' 'psoas' ...
    'quadfem' 'recfem' 'sart' 'semimem' 'semiten' 'soleus' 'tfl' 'tibant'...
    'tibpost' 'vasint' 'vaslat' 'vasmed'};
%
% %creates muscle amplitude files for trial
% muscfile=[musc_amp_path   trial '/' trial '-MUSC.inp'];
% if exist(muscfile,'file') ~= 2
%     disp('creating amplitute card')
%     MUSC=MUSC_AMP_CARDjw(path,patient,trial);
% end
trials ={ 'ngait_og1', 'stairdownf2', 'stairupf1', 'stepupf2', 'lungef1'}; %

for o=1:length(trials)
    trial=trials{o};
    force_file = [force_path trial '.mot'];
    kin_F_file = [ik_path trial '_F_ik.mot'];
    % muscfile=[musc_amp_path    trial '/' trial  '_step1-MUSC.inp'];
    % if exist(muscfile,'file') ~= 2
    %     disp('creating amplitute card')
    [ initial_time, final_time ] = gettrialtime(kin_F_file,force_file,patient, trial);
    time.(trial)  = getstoptimes(Path, patient,trial, initial_time, final_time);
    
    [TF, grfz] = MUSC_AMP_CARDjw(Path,patient,trial,initial_time,final_time,(time.(trial))); %each step
    
    totalforce.(trial)=TF;
    GRF.(trial)=grfz;
    % end
    % T = linspace(initial_time,final_time, size(MUSC,1));
    % plot(T,MUSC(:,[36:37 45:47]), 'LineWidth', 4)
    % legend(Musc_amp_list{[36:37 45:47]})
    
    
end
fluropath=[Path '/jw/FLURO/'];
save([musc_amp_path    '/totalforce.mat'],'totalforce');
save([musc_amp_path    '/GRF.mat'],'GRF');
save([fluropath 'steptimes.mat'],'time');

kneeflex.lungef1.step1=1.4347;
kneeflex.lungef1.step2=1.4927;
kneeflex.lungef1.step3=1.9069;
kneeflex.lungef1.step4=1.1549;

kneeflex.ngait_og1.step1=0.1481;
kneeflex.ngait_og1.step2=0.1708;
kneeflex.ngait_og1.step3=0.0525;
kneeflex.ngait_og1.step4=0.6113;

kneeflex.stairdownf2.step1=0.5095;
kneeflex.stairdownf2.step2=0.3622;
kneeflex.stairdownf2.step3=0.9053;
kneeflex.stairdownf2.step4=1.2064;

kneeflex.stairupf1.step1=1.0365;
kneeflex.stairupf1.step2=0.7210;
kneeflex.stairupf1.step3=0.2129;
kneeflex.stairupf1.step4=0.4262;

kneeflex.stepupf2.step1=1.0428;
kneeflex.stepupf2.step2=0.6037;
kneeflex.stepupf2.step3=0.5449;
kneeflex.stepupf2.step4=1.0318;

save([musc_amp_path    '/kneeflex.mat'],'kneeflex');

%hold conncetors constant. -ignore patella- doesnt work
%f=kin(:,[1:6 13:24]);%

%'R_PF_EXT_IE', 'R_PF_EXT_ML', 'R_PF_EXT_AP', 'R_PF_EXT_IE', 'R_PF_EXT_ML', 'R_PF_EXT_AP', ...






%%
load([musc_amp_path    'totalforce.mat']);
loads=struct;

% 
% loads.A.lungef1=[20000];
% loads.A.ngait_og1=[15000 15000 15000 8000] ;
% loads.A.stairdownf2=[10000 10000 20000 20000];
% loads.A.stairupf1=[15000 20000 15000 20000];
% loads.A.stepupf2=20000;
% 
% loads.P.lungef1=-20000;
% loads.P.ngait_og1=[-8000 -8000 -8000 -5000];
% loads.P.stairdownf2=[-3000 -3000 -8000 -8000];
% loads.P.stairupf1=[-8000 -8000 -8000 -4000];
% loads.P.stepupf2=-[15000 4000 4000 20000] ;
% 
% loads.I.lungef1=-[30000 30000 30000 30000 ];
% loads.I.ngait_og1=-[30000 30000 50000 50000 ];
% loads.I.stairdownf2=-[30000 30000 50000 30000 ];
% loads.I.stairupf1=-[50000 50000 30000 20000 ];
% loads.I.stepupf2=-[50000 30000 30000 30000 ];
% 
% loads.E.lungef1=[60000 50000 30000 60000 ];
% loads.E.ngait_og1=[30000 30000 50000 60000 ];
% loads.E.stairdownf2=[50000 30000 50000 30000 ];
% loads.E.stairupf1=[80000 80000 30000 30000 ];
% loads.E.stepupf2=[80000 30000 30000 80000 ];


loads.A.lungef1=[500];
loads.A.ngait_og1=500;% [15000 15000 15000 15000] ;
loads.A.stairdownf2=500;%[10000 10000 20000 20000];
loads.A.stairupf1=500;%[15000 20000 15000 20000];
loads.A.stepupf2=500;

loads.P.lungef1=-500;
loads.P.ngait_og1=-500;
loads.P.stairdownf2=-500;%[-5000 -5000 -10000 -10000];
loads.P.stairupf1=-500;%[-10000 -10000 -10000 -5000];
loads.P.stepupf2=-500;%[15000 4000 4000 20000] ;

loads.I.lungef1=-20000;%[30000 30000 30000 30000 ];
loads.I.ngait_og1=-20000;%-[30000 30000 50000 50000 ];
loads.I.stairdownf2=-20000;%-[30000 30000 50000 30000 ];
loads.I.stairupf1=-20000;%-[50000 50000 30000 20000 ];
loads.I.stepupf2=-20000;%-[50000 30000 30000 30000 ];

loads.E.lungef1=20000;%[60000 50000 30000 60000 ];
loads.E.ngait_og1=20000;%[30000 30000 50000 60000 ];
loads.E.stairdownf2=20000;%[50000 30000 50000 30000 ];
loads.E.stairupf1=20000;%[80000 80000 30000 30000 ];
loads.E.stepupf2=20000;%[80000 30000 30000 80000 ];

interval_time=0.025;
data={};

parameters=struct;
params=xlsread([Path patient '/INPUTS/parameters.xlsx'],1);
for i=1:7 %implant no
     impno=num2str(i);
    for j=1:4 %step in trial
        switch j
            case 1
                ind='A';
                case 2
                ind='E';
                case 3
                ind='I';
            case 4
                 ind='P';
        end
        
        for k=1:5 %trials
               switch k
            case 1
                tri='lungef1';
                case 2
                tri='ngait_og1';
                case 3
               tri='stairdownf2';
            case 4
                 tri='stairupf1';
                 case 5
                 tri='stepupf2';
               end
            C=i*5-4; %starting column (implant no)
            R=k*5-5+j; %starting row (trial)
           
          
        parameters.(['implant' impno ]).(tri).(ind)=params(R,[C:C+3]);
        
               
        end
    end
end
imp=[1:6];
%%%
for v=1:length(imp)
impno=imp(v);




implantno=num2str(impno);
a= fieldnames(parameters.implant3);
b=fieldnames(parameters.implant3.stepupf2);
% write the last field values in a single matrix
my_data = zeros(20,4);

for k=1:5 %trials
    for n=1:4 %step
        A=a{k};
        B=b{n};
        my_data(4*k-4+n,:) = parameters.implant5.(A).(B);
    end
end
% write the matrix to excell sheet

sf.imp1=1;
sf.imp2=1;
sf.imp3=1;
sf.imp4=1;
sf.imp5=1;

if exist([Path patient '/INPUTS/IMP' implantno '/FORTRAN' ], 'dir')==0
    mkdir([Path patient '/INPUTS/IMP' implantno '/FORTRAN' ])
end

%NON FORTRAN
for o=1:length(trials)
    trial=trials{o};
    %Stop input files
    names={'A','P','E','I','initial'};
    for step=1:4 %no steps
        num=num2str(step);
        for p=1:length(names)
            name=['step' num '_' names{p} ];
            
          
            
            file=[Path patient '/genericinput/genericIMP/genericjw_stepX_' names{p} '.inp'];
            fin_r = fopen(file,'r');
            r=1;
            while ~feof(fin_r) %Do until the end of the file
                temp = fgetl(fin_r);
                if contains(temp,'name')==1 || contains(temp,'trial')==1 || contains(temp,'GRFX')==1 ...
                        || contains(temp,'ALOAD')==1 || contains(temp,'PLOAD')==1 ||contains(temp,'GRFX')==1  || contains(temp,'MUSCX')==1   ...
                        || contains(temp,'ILOAD')==1 || contains(temp,'ELOAD')==1 ||  ...
                        contains(temp,'KINX')==1 || contains(temp,'MUSCX')==1 || contains(temp,'intTime')==1
                    temp = regexprep(temp, 'name', [patient trial]);
                    temp = regexprep(temp, 'trial', trial);
                    temp = regexprep(temp, 'KINX', ['KIN' num]);
                    temp = regexprep(temp, 'GRFX', ['GRF' num]);
                    temp = regexprep(temp, 'MUSCX', ['MUSC' num]);
                    if length(loads.A.(trial))==1
                        A=num2str(loads.A.(trial));
                    else
                        A=num2str(loads.A.(trial)(step));
                    end
                    temp = regexprep(temp, 'ALOAD', A);
                    
                    if length(loads.P.(trial))==1
                        P=num2str(loads.P.(trial));
                    else
                        P=num2str(loads.P.(trial)(step));
                    end
                    temp = regexprep(temp, 'PLOAD', P);
                    
                    if length(loads.E.(trial))==1
                        E=num2str(loads.E.(trial));
                    else
                        E=num2str(loads.E.(trial)(step));
                    end
                    temp = regexprep(temp, 'ELOAD', E);
                    
                             if length(loads.I.(trial))==1
                        I=num2str(loads.I.(trial));
                    else
                        I=num2str(loads.I.(trial)(step));
                    end
                    temp = regexprep(temp, 'ILOAD', I);
                    
          
                    temp = regexprep(temp, 'intTime', num2str(interval_time));
                    
                    data{r}=temp;
                    r=r+1;
                    %
                else
                    data{r}=temp;
                    r=r+1;
                end
            end
            
            fin_w1 = fopen([Path patient '/INPUTS/IMP_base/' trial name '.inp'],'w+');
            for i=1:r-1
                fprintf(fin_w1,'%s',data{i});
                fprintf(fin_w1,'\n');
            end
            fclose('all');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%Fortran files%%%%
        load([Path patient '/genericinput/Initialdisp/DISP_ROT_IMP' implantno '.mat']);
        
        names2={'A','P','E','I', 'initial'};
        
        musc.initial={   'gaslat' 'gasmed' ...
                'recfem'  'vasint' 'vaslat' 'vasmed' ...
                'bflh' 'bfsh' 'grac'   'semimem' 'semiten' 'tfl' 'sart'};
            
            musc.A={  'gaslat' 'gasmed' ...
                'bflh' 'bfsh' 'grac'   'semimem' 'semiten' 'tfl' 'sart'};
            
            musc.P=    {    'recfem'  'vasint' 'vaslat' 'vasmed'};
            
                      musc.I=  {  'gaslat' ...
                'recfem'  'vasint' 'vaslat' 'vasmed' ...
                'bflh' 'bfsh'    'tfl' }; %External rotators lateral hammy
            
            musc.E={   'gasmed' ...
                'recfem'  'vasint' 'vaslat' 'vasmed' ...
                'grac'   'semimem' 'semiten'  'sart'}; %internal rotators =medial

            
        for p=1:length(names2)
            name=['step' num '_' names2{p} ];
            
            
            file=[Path patient '/genericinput/fortran/genericjw_stepX_' names2{p} '.inp'];
            if impno==6 || impno==7
                    file=[Path patient '/genericinput/fortran/genericjw_stepX_' names2{p} '_patient.inp'];
            end
                
              if strcmp(trial,'lungef1')==1 %|| strcmp(trial,'ngait_og1')==1 
                if strcmp(name,'step3_initial')==1
                    file= insertAfter(file,'initial','newload');
                end
              end
            
            if p==5
                file2=[Path patient '/genericinput/fortran/genericfortran-initial.f'];
                

                
            elseif p==1 || p==2
                file2=[Path patient '/genericinput/fortran/genericfortranAP.f'];
            else
                file2=[Path patient '/genericinput/fortran/genericfortranIE.f'];
                
            end
            
            fin_r = fopen(file,'r');
            fin_r2 = fopen(file2,'r');
            
            interval_time=0.025;%TotTime/100;
            r=1;
            
            while ~feof(fin_r) %Do until the end of the file
                temp = fgetl(fin_r);
                if contains(temp,'name')==1 || contains(temp,'trial')==1 || contains(temp,'GRFX')==1 ...
                        || contains(temp,'ALOAD')==1 || contains(temp,'PLOAD')==1 ||contains(temp,'GRFX')==1  || contains(temp,'MUSCX')==1   ...
                        || contains(temp,'ILOAD')==1 || contains(temp,'ELOAD')==1 || contains(temp,'frac_')==1 ||...
                        contains(temp,'TotalF')==1 || contains(temp,'MUSCX')==1 || contains(temp,'intTime')==1 ...
                        || contains(temp,'InsertX')==1
                    temp = regexprep(temp, 'name', [patient trial]);
                    temp = regexprep(temp, 'trial', trial);
                    temp = regexprep(temp, 'KINX', ['KIN' num]);
                    temp = regexprep(temp, 'GRFX', ['GRF' num]);
                    temp = regexprep(temp, 'MUSCX', ['MUSC' num]);
                    temp = regexprep(temp, 'InsertX', ['Insert' implantno]);
                    if length(loads.A.(trial))==1
                        A=(loads.A.(trial));
                    else
                        A=(loads.A.(trial)(step));
                    end
                    temp = regexprep(temp, 'ALOAD', num2str(A));
                    
                    if length(loads.P.(trial))==1
                        P=(loads.P.(trial));
                    else
                        P=(loads.P.(trial)(step));
                    end
                    temp = regexprep(temp, 'PLOAD', num2str(P));
                    
                    if length(loads.E.(trial))==1
                        E=(loads.E.(trial));
                    else
                        E=(loads.E.(trial)(step));
                    end
                    temp = regexprep(temp, 'ELOAD', num2str(E));
                    %
                          if length(loads.I.(trial))==1
                        I=(loads.I.(trial));
                    else
                        I=(loads.I.(trial)(step));
                    end
                    temp = regexprep(temp, 'ILOAD', num2str(I));
                  
           
                    temp = regexprep(temp, 'intTime', num2str(interval_time));
                    temp = regexprep(temp, 'TotalF', num2str(totalforce.(trial).(['step' num]).(names2{p})));
                    
                    list= musc.(names{p});
                    
                    for m=1:length(list)
                        fract=totalforce.(trial).fract.(['step' num]).(names{p}).(list{m});
                        temp = regexprep(temp, ['frac_' list{m}] , num2str(-1*fract));
                    end
                    
                    data{r}=temp;
                    r=r+1;
                    %
                else
                    data{r}=temp;
                    r=r+1;
                end
            end
            
            if exist([Path patient '/INPUTS/FORTRAN/IMP' implantno  ], 'dir')==0
                mkdir([Path patient '/INPUTS/FORTRAN/IMP' implantno ])
            end
            
            fin_w1 = fopen([Path patient '/INPUTS/FORTRAN/IMP' implantno '/' trial name '_IMP' implantno '-FT.inp'],'w+');
            for i=1:r-1
                fprintf(fin_w1,'%s',data{i});
                fprintf(fin_w1,'\n');
            end
            
            r=1;
            
            
            %%F file
            %%%%%%%%%%%%%%%%%%%%%
            while ~feof(fin_r2) %Do until the end of the file
                temp = fgetl(fin_r2);
                if contains(temp,'TotalF')==1 || contains(temp,'APLOAD')==1 ||   ...
                        contains(temp,'IELOAD')==1 || contains(temp,'TRIALNAMESTEPNO')==1 ||...
                        contains(temp,'Ki_AP')==1 || contains(temp,'Kp_AP')==1 || ...
                        contains(temp,'Ki_IE')==1 || contains(temp,'Kp_IE')==1 || ...
                        contains(temp,'initialrot')==1  || contains(temp,'IMPno')==1
                    temp = regexprep(temp, 'TotalF', num2str(totalforce.(trial).(['step' num]).(names2{p})));
                    temp = regexprep(temp, 'IMPno', ['IMP' num2str(implantno)]);
                    try
                        temp = regexprep(temp, 'initialrot', num2str(disp_rot.IE.([trial 'step' num '_initial_IMP' implantno]))); %should be same for all AIEP
                    catch
                    end
                    
                    
                    
                    if contains(names{p},'initial')==0
                        ind=parameters.(['implant' implantno]).(trial).(names{p})(step);
                        
                        switch   ind
                            case 0
                                Kp=10 ; Ki=50;
                            case 1
                                Kp=.05;   Ki=.500;
                            case 2
                                Kp=2;   Ki=100;
                            case 3
                                Kp=.02;   Ki=.2;
                                 case 3.2
                                Kp=.02;   Ki=20;
                                 case 3.1
                                Kp=.02;   Ki=10;
                            case 3.5
                                Kp=.02;   Ki=10;
                                case 3.9
                                Kp=.02;   Ki=30;
                            case 4
                                Kp=.01;   Ki=1;
                                 case 4.1
                                Kp=.01;   Ki=100;
                            case 5
                                Kp=.01;   Ki=10;
                                 case 5.5
                                Kp=.01;   Ki=30;
                            case 6
                                Kp=.025;   Ki=5;
                            case 7
                                Kp=.05;   Ki=5;
                            case 7.5
                                Kp=.05;   Ki=10;
                                    case 7.8
                                Kp=.04;   Ki=20;
                                   case 7.9
                                Kp=.05;   Ki=50;
                            case 8
                                Kp=.075;   Ki=.75;
                            case 9
                                Kp=.1;   Ki=1;
                            case 9.5
                                Kp=.1;   Ki=5;
                            case 10
                                Kp=.1;   Ki=15;
                                 case 10.5
                                Kp=.1;   Ki=50;
                            case 11
                                Kp=.5;   Ki=5;
                                   case 11.5
                                Kp=.5;   Ki=20;
                            case 12
                                Kp=.5;   Ki=50;
                            case 13
                                Kp=1;   Ki=1;
                                case 13.2
                                Kp=1;   Ki=5;
                            case 13.5
                                Kp=1;   Ki=20;
                            case 14
                                Kp=1;   Ki=50;
                            case 15
                                Kp=2;   Ki=10;
                                 case 15.5
                                Kp=2;   Ki=20;
                            case 16
                                Kp=2;   Ki=100;
                                 case 16.5
                                Kp=3;   Ki=15;
                            case 17
                                Kp=3;   Ki=200;
                                  case 17.1
                                Kp=3;   Ki=9;
                            case 17.5
                                Kp=3;   Ki=50;
                            case 18
                                Kp=4;   Ki=100;
                            case 18.5
                                Kp=4;   Ki=12;
%                                   case 18.5
%                                 Kp=5;   Ki=1;
                          
                            case 19
                                Kp=5;   Ki=250;
                            case 19.5
                                Kp=5;   Ki=15;
                            case 20
                                Kp=10;   Ki=10;
                                   case 20.2
                                Kp=10;   Ki=50;
                            case 20.5
                                Kp=10;   Ki=500;
                            case 21
                                Kp=7;   Ki=300;
                            case 22
                                Kp=8;   Ki=200;
                                     case 22.5
                                Kp=7.5;   Ki=350;
                            case 23
                                Kp=.5;   Ki=500;
                                  case 23.2
                                Kp=.25;   Ki=100;
                                   case 23.5
                                Kp=.5;   Ki=200;
                                  case 24
                                Kp=20;   Ki=20;
                        end
                        
                    else
                        Kp=0;
                        Ki=0;
                        
                    end
                    
                    if impno==7
                    Kp=0;
                    Ki=0;
                    end
                    
                    if p==1
                        
                        temp = regexprep(temp, 'APLOAD', num2str(A));
                        temp = regexprep(temp, 'IELOAD', '0.0');
                    elseif p==2
                        temp = regexprep(temp, 'APLOAD', num2str(abs(P)));
                        temp = regexprep(temp, 'IELOAD', '0.0');
                    elseif p==3
                        temp = regexprep(temp, 'APLOAD', '0.0');
                        temp = regexprep(temp, 'IELOAD', num2str(E));
                    elseif p==4
                        temp = regexprep(temp, 'APLOAD', '0.0');
                        temp = regexprep(temp, 'IELOAD', num2str(abs(I)));
                    end
                    %         temp = regexprep(temp, 'ALOAD', num2str(loads.A.(trial)));
                    %         temp = regexprep(temp, 'PLOAD', num2str(loads.P.(trial)));
                    %         temp = regexprep(temp, 'ILOAD', num2str(loads.I.(trial)));
                    %         temp = regexprep(temp, 'ELOAD', num2str(loads.E.(trial)));
                    temp = regexprep(temp, 'TRIALNAMESTEPNO', [trial 'step' num]);
                    temp = regexprep(temp, 'TYPE', ['_' names2{p}]);
                    
                    temp = regexprep(temp, 'Kp_AP', num2str(Kp));
                    temp = regexprep(temp, 'Ki_AP', num2str(Ki));
                    temp = regexprep(temp, 'Kp_IE', num2str(Kp));
                    temp = regexprep(temp, 'Ki_IE', num2str(Ki));
                    data{r}=temp;
                    r=r+1;
                    %
                else
                    data{r}=temp;
                    r=r+1;
                end
            end
            
            
            fin_w2 = fopen([Path patient '/INPUTS/FORTRAN/IMP' implantno '/' 'vuamp-' trial name '_IMP' implantno '-FT.f'],'w+');
            for i=1:r-1
                fprintf(fin_w2,'%s',data{i});
                fprintf(fin_w2,'\n');
            end
            
            %
            % if exist ([path patient '/RESULTS/' trial '_' name '/']) ==0;
            %     mkdir([path patient '/RESULTS/' trial '_' name '/']);
            %
            %
            % fin_w3 = fopen([path patient '/RESULTS/' trial '_' name '/' 'diag-ap.txt'],'w+');
            % fin_w4 = fopen([path patient '/RESULTS/' trial '_' name '/' 'diag-ie.txt'],'w+');
            % fin_w5 = fopen([path patient '/RESULTS/' trial '_' name '/' 'target_AP.inp'],'w+');
            % fin_w6 = fopen([path patient '/RESULTS/' trial '_' name '/' 'target_IE.inp'],'w+');
            %
            % end
            
            fclose('all');
        end
        
        
    end
    fclose('all');
end
end



%%
disp('test7 -%Creates GS Kin Axes file')
%Creates GS Kin Axes file

axes_list = {'R_TF_EXT_IE', 'R_TF_EXT_ML', 'R_TF_EXT_AP',  ...
    'R_PF_EXT_IE', 'R_PF_EXT_ML', 'R_PF_EXT_AP', ...
    'R_Hip_IE', 'R_Hip_ML', 'R_Hip_AP',  ...
    'R_Ank_IE', 'R_Ank_ML', 'R_Ank_AP'};
joint_list = {'R_TF_EXT', 'R_PF_EXT', 'R_Hip',  'R_Ank'};

%Reference Frame Nodes

%node numbers for reference nodes.
for i = 1:5
    switch i
        case 1; I= 100;
        case 2; I=200;
        case 3;  I=300;
        case 4; I=800;
        case 5;  I=900;
    end
    for j=1:4
        N(4*i+j-4)=[I+j-1]; %100-103 etc
    end
end
fout1=fopen([Path patient '/LtModel/GS_Kin_Axes_' name '.inp'],'w+');
% fprintf(fout1,'**\n*NODE, nset=ref_frames\n');
% % no=n;
% % no(1:4:end,:) = []; %dont want origin node already defined in frame input file
% for i = 1:20 %defines origin?ref frame nodes nodes
%     fprintf(fout1,'%1.0f, %1.5f, %1.5f, %1.5f\n',[N(i) n(i,2:4)]);
% end
fprintf(fout1,'**\n*mpc');
for i = 2:4
    fprintf(fout1,'\nbeam, '); %beams orgin reference nodes to each other for each bone
    fprintf(fout1, '%1.0f, %1.0f', [N(i) 100]);
end
for i = 6:8
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) 200]);
end
for i = 10:12
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) 300]);
end
for i = 14:16
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) 800]);
end
for i = 18:20
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) 900]);
end
N=[104 105 106 107];
for i = 1:4
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) 100]);
end
N=[204 205 206 207];
for i = 1:4
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) 200]);
    
end

for i = 1:4
    fprintf(fout1,['\n**\n**\n**Intersection nodes to define floating AP axes - ' joint_list{i} '\n*NODE\n']);
    switch joint_list{i} %prints node no.s and coordinates
        %Cor1=top bone %Cord2=btm bone node
        case 'R_TF_EXT'; Node = [150 151]; Cor1 = inter_rtf1(1,:); Cor2 = inter_rtf2(1,:);
        case 'R_PF_EXT'; Node = [250 251]; Cor1 = inter_rpf1(1,:); Cor2 = inter_rpf2(1,:);
        case 'R_Hip'; Node = [350 351]; Cor1 = inter_rhip1(1,:); Cor2 = inter_rhip2(1,:);
        case 'R_Ank'; Node = [450 451]; Cor1 = inter_rank1(1,:); Cor2 = inter_rank2(1,:);
            % Cor1(1)=Cor1(1)*-1; Cor2(1)=Cor2(1)*-1;
    end
    fprintf(fout1, '%1.0f, %1.5f, %1.5f, %1.5f', [Node(1) Cor1]);
    fprintf(fout1, '\n%1.0f, %1.5f, %1.5f, %1.5f', [Node(2) Cor2]);
end


for i = 1:4
    %     switch joint_list{i}
    %         case 'R_TF_EXT'; Node = [150 151]; Element = [150 151 152];  Ori = [200 100]; %%element =axes
    %         case 'R_PF_EXT'; Node = [250 251]; Element = [250 251 252]; Ori = [300 100];
    %         case 'R_Hip'; Node = [350 351]; Element = [350 351 352]; Ori = [800 4]; %4=GS node femur hip/superior
    %         case 'R_Ank'; Node = [450 451]; Element = [450 451 452]; Ori = [7 900]; %7=GS tibiba  inferior node
    %     end
    
    switch joint_list{i}
        case 'R_TF_EXT'; Node = [150 151]; Element = [150 151 152];  Ori = [200 100]; %%element =axes
        case 'R_PF_EXT'; Node = [250 251]; Element = [250 251 252]; Ori = [    300 100 ];
        case 'R_Hip'; Node = [350 351]; Element = [350 351 352]; Ori = [    4 800]; %4=top GS node femur hip/superior
        case 'R_Ank'; Node = [450 451]; Element = [450 451 452]; Ori = [   7 900]; %7=GS tibiba  inferior node
    end
    
    fprintf(fout1,['\n**\n**\n********** ' joint_list{i} ' Kinematic Axes\n**\n*ELEMENT,TYPE=CONN3D2,ELSET=axis_' joint_list{i} '_IE\n']);
    fprintf(fout1, '%1.0f, %1.0f, %1.0f', [Element(1)   Node(1) Ori(1) ]); %IE axis =1st element between top origin GS node and top origin AP1
    fprintf(fout1,['\n*ELEMENT,TYPE=CONN3D2,ELSET=axis_' joint_list{i} '_ML\n']);
    fprintf(fout1, '%1.0f, %1.0f, %1.0f', [Element(3)  Node(2) Ori(2)   ]);
    fprintf(fout1,['\n*ELEMENT,TYPE=CONN3D2,ELSET=axis_' joint_list{i} '_AP\n']);
    fprintf(fout1, '%1.0f, %1.0f, %1.0f', [Element(2)  Node(2) Node(1)]); %flipped all
    
    
    %%original
    %IE ML AP
    %     if i == 1  %TF Joints
    %         for j = 1:3
    %             k = 3*(i-1)+j;
    %             fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
    %             if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(1) Ori(1)+1 Ori(1)]); %[150 202 200]
    %             elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(2) Ori(2)+2 Node(2)]);
    %             elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(2) Ori(2) Node(1)]);
    %             end
    %         end
    %     elseif i == 2  %PF Ori = [  300 100];
    %         for j = 1:3
    %             k = 3*(i-1)+j;
    %             fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
    %             if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Ori(1)+1 Ori(1)  Node(1)]); %Ori(1)+1 Ori(1)
    %             elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Node(2) Ori(2)+2 Ori(2)]); %
    %             elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Ori(2) Node(1) Node(2)]);
    %             end
    %         end
    %     elseif i == 3  %Hip Joints
    %         for j = 1:3
    %             k = 3*(i-1)+j;
    %             fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
    %             if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(1) Ori(1)+1 Node(1)]);
    %             elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(2) Ori(2)+2 Node(2)]);
    %             elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(2) Ori(2) Node(1)]);
    %             end
    %         end
    %     elseif i == 4  %Ankle Joints
    %         for j = 1:3
    %             k = 3*(i-1)+j;
    %             fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
    %             if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(1) Ori(1)+1 Node(1)]);
    %             elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(2) Ori(2)+2 Node(2)]);
    %             elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(1) Ori(2) Node(2)]);
    %             end
    %         end
    %     end
    % end
    
    %     orientation of axes
    %coordinates of the three points a, b, and, optionally, c (the origin)
    %cylindrical coordinate system with a and b on polar axis and c at origin
    %IE =SI translation, int/ext rotation pointing up +vs=int rotation
    %ML=+ve to the right, Flexion =+ve
    %AP= +ve is posterior and adduction
    %working
    
    
    %
    if i == 1 || i==2 || i==4 %|| i==3 %TF Joints
        for j = 1:3
            k = 3*(i-1)+j;
            fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
            if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Ori(1) Ori(1)+1 Node(1)]); %tried IE FLIPPED[150 202 200]  orig [Ori(1)+2 Ori(1)  Node(1)])
            elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Node(2) Ori(2)+2  Ori(2)   ]); %ML Node(2) Ori(2)+3 Ori(2)
            elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f',  [  Ori(2) Node(2)   Node(1)]); %AP tried fliiping didnt work.  [  Ori(2) Node(2)   Node(1)])
            end
        end
        %WORKING PATELLA I THINK %IE ML AP
    elseif i == 2  %PF Ori = [  300 100];
        for j = 1:3
            k = 3*(i-1)+j;
            fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
            if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Node(1)  Ori(1)+1  Ori(1) ]); %Node(1)   Ori(1) Ori(1)+1 %change to +2? no because
                %diff coord system??
            elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [  Node(2)  Ori(2)+2  Ori(2)]);  %changed to +3  Node(2) Ori(2)+3 Ori(2)
            elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Node(1) Node(2) Ori(2)]); % SWITCHED FROM
            end
        end
        
        %     elseif i == 3  %Hip Joints Ori = [ 104 800 ];
        %         for j = 1:3 %IE ML AP
        %             k = 3*(i-1)+j;
        %             fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
        %             if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Ori(1)+2   Ori(1) Node(1)]); %IE +y %femur 150 201 200
        %             elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [   Node(2) Ori(2)+3  Ori(2) ]); %flip since decreasing angle =flexion ML 100 202 151
        %             elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Ori(2) Node(1)   Node(2)]); %AP
        %             end
        %         end
        
        %best so far
    elseif i == 3  %Hip Joints
        for j = 1:3
            k = 3*(i-1)+j;
            fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
            if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(1) Ori(1)+1 Node(1)  ]); %IE
            elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Node(2) Ori(2)+2 Ori(2) ]); %ML
            elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(2) Ori(2) Node(1)]); %AP
                %                    if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(1) Ori(1)+1 Node(1) ]); %IE
                %             elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Node(2) Ori(2)+2 Ori(2)  ]);
                %             elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Node(2) Ori(2) Node(1)]); %Node(2) Ori(2) Node(1)
            end
        end
        %
        
    elseif i == 4  %Ankle Joints
        for j = 1:3
            k = 3*(i-1)+j;
            fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
            if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Ori(1) Ori(1)+1 Node(1)]); %IE tibia
            elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Node(2) Ori(2)+2  Ori(2)   ]); %ML talus
            elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [ Ori(2) Ori(2)+2  Node(1)]); %AP
            end
        end
    end
    %
end


Node = [150 151 250 251  350 351 450 451 ];

fprintf(fout1,'\n**\n**\n**\n**\n********\n**\n*ELEMENT, TYPE=MASS, ELSET=CONN_KIN_MASS');
for i = 1:8
    fprintf(fout1, '\n%1.0f, %1.0f', [i Node(i)]);
end

fprintf(fout1,'\n*ELEMENT, TYPE=ROTARYI, ELSET=CONN_KIN_ROT');
for i = 1:8
    fprintf(fout1, '\n%1.0f, %1.0f', [i+8 Node(i)]);
end

fprintf(fout1,'\n**\n*MASS, ELSET=CONN_KIN_MASS\n1e-5\n*ROTARY INERTIA, ELSET=CONN_KIN_ROT\n0.001, 0.001, 0.001\n**');
fclose(fout1);
%%
disp('test18')
%writes abaqus input file and output file replacing patient name, final time and output
%time increments

file=[Path patient '/genericinput/genericjw.inp'];
fin_r = fopen(file,'r');
TotTime=(final_time-initial_time);%from SLOW amp cards

interval_time=0.02;%TotTime/100;
r=1;

while ~feof(fin_r) %Do until the end of the file
    temp = fgetl(fin_r);
    if contains(temp,'name')==1
        temp = regexprep(temp, 'name', [patient trial]);
        data{r}=temp;
        r=r+1;
    elseif contains(temp,'TotTime')==1
        temp = regexprep(temp, 'TotTime', num2str(TotTime));
        data{r}=temp;
        r=r+1;
    elseif contains(temp,'intTime')==1
        temp = regexprep(temp, 'intTime', num2str(interval_time));
        data{r}=temp;
        r=r+1;
    elseif contains(temp,'trial')==1
        temp = regexprep(temp, 'trial', trial);
        data{r}=temp;
        r=r+1;
    else
        data{r}=temp;
        r=r+1;
    end
end


fin_w1 = fopen([Path patient '\' name '.inp'],'w+');
for i=1:r-1
    fprintf(fin_w1,'%s',data{i});
    fprintf(fin_w1,'\n');
end

fclose(fin_w1);

fclose('all');

%origianl
% loads.A.lungef1=2*[7500];
% loads.A.ngait_og1=2*[7500 7500 7500 3000] ;
% loads.A.stairupf1=2*[7500];
% loads.A.stairdownf2=2*7500;
% loads.A.stepupf2=2*7500;
% 
% loads.P.lungef1=2*-7500;
% loads.P.ngait_og1=2*[-3000 -3000 -3000 -1000];
% loads.P.stairupf1=2*[-3000 -3000 -3000 -1000];
% loads.P.stairdownf2=2*[-1500 -1500 -3000 -1500];
% loads.P.stepupf2=-2*[5000 1500 1500 7500] ;
% 
% loads.I.lungef1=-30000;
% loads.I.ngait_og1=-30000;
% loads.I.stairupf1=-30000;
% loads.I.stairdownf2=-30000;
% loads.I.stepupf2=-30000;
% 
% loads.E.lungef1=[30000];
% loads.E.ngait_og1=[30000 ];
% loads.E.stairupf1=[30000 ];
% loads.E.stairdownf2=[30000 ];
% loads.E.stepupf2=[30000];



%  
% parameters.implant1.lungef1.A=   [17 1 21 21 ]; %trying .5 1
% parameters.implant1.lungef1.E=   [7.5 0 15 19.5];
% parameters.implant1.lungef1.I=   [11 0 6 15 ];
% parameters.implant1.lungef1.P=   [7 1 7 6];
%  
% parameters.implant1.ngait_og1.A=   [21 21 23 20.5 ]; %try 20
% parameters.implant1.ngait_og1.E=   [17.5 17.5 17.5 7.5 ];
% parameters.implant1.ngait_og1.I=  [17.5 15 14 10 ];
% parameters.implant1.ngait_og1.P=  [20.5 17 21 9.5];
%  
% parameters.implant1.stairdownf2.A=  [20.5 16 17 21 ]; %try 19
% parameters.implant1.stairdownf2.E= [16 15 14 15 ];
% parameters.implant1.stairdownf2.I= [15 9 17.5 15 ];
% parameters.implant1.stairdownf2.P=  [7.5 10 21 10];
%  
% parameters.implant1.stairupf1.A= [17 20.5 21 20.5];
% parameters.implant1.stairupf1.E= [13 15 15 16 ];
% parameters.implant1.stairupf1.I= [10 10 7.5 15 ];
% parameters.implant1.stairupf1.P= [3 4 21 7.9 ];
%  
% parameters.implant1.stepupf2.A= [21 20.5 20.5 21 ];
% parameters.implant1.stepupf2.E= [15 10 16 15 ];
% parameters.implant1.stepupf2.I= [10 17.5 18 10 ];
% parameters.implant1.stepupf2.P= [3.5 5 23.5 5 ];
% 
% 
% %%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%
% parameters.implant2.lungef1.A=   [17 0 17 17  ];
% parameters.implant2.lungef1.E=   [7 0 11 19.5 ];
% parameters.implant2.lungef1.I=   [7 0 7 14];
% parameters.implant2.lungef1.P=   [7 0 7 6];
%  
% parameters.implant2.ngait_og1.A=   [21 21 23 20.5 ];
% parameters.implant2.ngait_og1.E=   [14 17 17.5 18.5 ];
% parameters.implant2.ngait_og1.I=  [14 14 14 14 ];
% parameters.implant2.ngait_og1.P=  [21 17 21 10 ];
%  
% parameters.implant2.stairdownf2.A= [20.5 23 19.5 22 ];
% parameters.implant2.stairdownf2.E=  [14 15 15.5 18.5 ];
% parameters.implant2.stairdownf2.I=  [14 15 10 14 ];
% parameters.implant2.stairdownf2.P=  [14 13 18.5 14 ];
%  
% parameters.implant2.stairupf1.A= [17 20.5 20.5 20.5 ];
% parameters.implant2.stairupf1.E= [13 18.5 14 13.5];
% parameters.implant2.stairupf1.I= [10 14 14 17.5 ];
% parameters.implant2.stairupf1.P= [3 3.5 21 7.9];
%  
% parameters.implant2.stepupf2.A= [21 20.5 20.5 21 ];
% parameters.implant2.stepupf2.E= [15 4.1 18 15];
% parameters.implant2.stepupf2.I= [14 17.5 17.5 13.5 ];
% parameters.implant2.stepupf2.P= [3.5 5 23.5 14 ];
% 
% 
% %%%%%%%%%%%%%%%%
% parameters.implant3.lungef1.A=   [17 1 17 21];
% parameters.implant3.lungef1.E=   [7.5 0 15 5];
% parameters.implant3.lungef1.I=   [11 0 11 15 ];
% parameters.implant3.lungef1.P=   [7 1 7 6];
%  
% parameters.implant3.ngait_og1.A=   [21  21  23  20.5 ];
% parameters.implant3.ngait_og1.E=   [17.5 17.5 16 15 ];
% parameters.implant3.ngait_og1.I=  [17.5 15 17.5 13.5 ];
% parameters.implant3.ngait_og1.P=  [21 17 21 11 ];
%  
% parameters.implant3.stairdownf2.A=  [20.5 23 12 22 ];
% parameters.implant3.stairdownf2.E= [17.5 15 13.5 15 ];
% parameters.implant3.stairdownf2.I= [15 9 17.5 15 ];
% parameters.implant3.stairdownf2.P=  [18.5 13 21 7];
%  
% parameters.implant3.stairupf1.A= [17 20.5 21 21];
% parameters.implant3.stairupf1.E= [13 15  15 17.5 ];
% parameters.implant3.stairupf1.I= [10 13.5 17.5 15 ];
% parameters.implant3.stairupf1.P= [3 4 22.5 17.5 ];
%  
% parameters.implant3.stepupf2.A= [21 20.5 20.5 21 ];
% parameters.implant3.stepupf2.E= [15 10 16 15];
% parameters.implant3.stepupf2.I= [10 17.5 18 10 ];
% parameters.implant3.stepupf2.P= [3.5 3 23.5 5 ];
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% parameters.implant4.lungef1.A=   [16 0 16 16  ];
% parameters.implant4.lungef1.E=   [7.5 0 9 7 ];
% parameters.implant4.lungef1.I=   [8 0 8 15];
% parameters.implant4.lungef1.P=   [7 0 7 5];
%  
% parameters.implant4.ngait_og1.A=   [19 19 18.5 21 ];
% parameters.implant4.ngait_og1.E=   [10 10 10 15 ];
% parameters.implant4.ngait_og1.I=  [10 10 10 3.5 ];
% parameters.implant4.ngait_og1.P= [21 21 21 13.5 ];
%  
% parameters.implant4.stairdownf2.A=  [21 20.5 21 21];
% parameters.implant4.stairdownf2.E=  [13.5 10 13.5 10 ];
% parameters.implant4.stairdownf2.I=  [10 9 10 10 ];
% parameters.implant4.stairdownf2.P=  [7 13 8 7];
%  
% parameters.implant4.stairupf1.A= [17 19 21 18.5 ];
% parameters.implant4.stairupf1.E= [15 15 10 10 ];
% parameters.implant4.stairupf1.I= [10 13.5 10 15 ];
% parameters.implant4.stairupf1.P= [3 4 19.5 18.5];
%  
% parameters.implant4.stepupf2.A= [17 20.5 20.5 21 ];
% parameters.implant4.stepupf2.E= [15 10 7 7.5 ];
% parameters.implant4.stepupf2.I= [10 17.5 15.5 10 ];
% parameters.implant4.stepupf2.P= [8 7.9 7.9 4 ];
%   
% parameters.implant5.lungef1.A=   [16 0 16 16  ];
% parameters.implant5.lungef1.E=   [7.5 0 9 7 ];
% parameters.implant5.lungef1.I=   [8 0 8 15];
% parameters.implant5.lungef1.P=   [7 0 7 5];
%  
% parameters.implant5.ngait_og1.A=   [19 19 18.5 21 ];
% parameters.implant5.ngait_og1.E=   [10 10 10 15 ];
% parameters.implant5.ngait_og1.I=  [10 10 10 3.5 ];
% parameters.implant5.ngait_og1.P= [21 21 21 13.5 ];
%  
% parameters.implant5.stairdownf2.A=  [21 20.5 21 21];
% parameters.implant5.stairdownf2.E=  [13.5 10 13.5 10 ];
% parameters.implant5.stairdownf2.I=  [10 9 10 10 ];
% parameters.implant5.stairdownf2.P=  [7 13 8 7];
%  
% parameters.implant5.stairupf1.A= [17 19 21 18.5 ];
% parameters.implant5.stairupf1.E= [15 15 10 10 ];
% parameters.implant5.stairupf1.I= [10 13.5 10 15 ];
% parameters.implant5.stairupf1.P= [3 4 19.5 18.5];
%  
% parameters.implant5.stepupf2.A= [17 20.5 20.5 21 ];
% parameters.implant5.stepupf2.E= [15 10 7 7.5 ];
% parameters.implant5.stepupf2.I= [10 17.5 15.5 10 ];
% parameters.implant5.stepupf2.P= [8 7.9 7.9 4 ];
