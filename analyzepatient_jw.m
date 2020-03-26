function []=analyzepatient_jw(path,patient,trial)
%%

patient = 'jw';
%
trial ='ngait_og2';
%path='C:/Users/Public/Google_Drive/Thesis/ACL_Pilot_data/Processing/';
path = ['C:/Users/Gracemcconnochie/Google_Drive/Thesis/ACL_Pilot_data/Processing/'];
cd(path)
name = [patient trial]; %used as name in analysis transform, make sure that works
%
import org.opensim.modeling.* % % Import OpenSim Libraries
version = 'V6'; %Changes name of HM GS Kin Axis File
threshold_m = 0; %Minimum threshold of muscle forces
threshmax = 1500; %Maximum threshold of muscle forces
musc_multi = 1; %Multiplier of muscle forces
MPFL = '-8';
LPFL = '-5';
vs = 7843; %Call number for HM Tibia superior node
vinf = 8225; %Call number for HM Tibia inferior node

trans_path = [path patient '/T-Matrix/']; % fill with path
scale_path = [path patient '/Scaling/'];
trc_path= [path patient '/MarkerData/'];
ik_path =[path patient '/IKResults/'];
state_path =[path patient '/StatesResults/'];
amp_path= [path patient '/Amplitudes/'];
force_path= [path patient '/ForceData/'];
ct_path = [path patient '/LtModel/'];
STO_path=[path patient '/STOResults/' ];
% Model
model_file = [scale_path patient '_Model.osim'];
model = Model(model_file);



%%
disp('test0');
%creates necessary setup files

force_file = [force_path trial '.mot'];
kin_F_file = [ik_path trial '_F_ik.mot']; % fill with output file of IK (.mot)
if exist(kin_F_file,'file') ~= 2
    disp('Running IK full cycle analysis')
    setupAndRunIKSingleCycle( patient,path, trial)
end
[ initial_time, final_time ] = gettrialtime(kin_F_file,force_file,patient, trial);

%run cropped inverse kinematics
kin_file = [ik_path trial '_ik.mot']; % fill with output file of IK (.mot)
if exist(kin_file,'file') ~= 2
    disp('Running IK cropped cycle analysis')
    setupAndRunIKSingleCycle( patient,path, trial,initial_time, final_time )
end

state_file = [state_path trial '_StatesReporter_states.sto'];
% output of the StateReporter Analyze tool
if exist(state_file,'file') ~= 2
    disp('Running States analysis')
    setupAndRunStatesSingleCycle( patient,path, trial, initial_time, final_time)
end

STOfile=[STO_path trial '_StaticOptimization_force.sto'];
if exist(STOfile,'file') ~= 2
    disp('Running STO analysis')
    setupAndRunSTOSingleCycle(patient, path, trial, initial_time, final_time)
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
    load([trans_path name '.mat'])
    disp('loaded transformation matrices')
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
    
    
    save([trans_path name '.mat'],'TM') %'-struct',
    
   %load([trans_path name '.mat'],'T*');
end


disp('test2')


%%
disp('test3')
%Bone Coordinate Transformation
%Creates Bone inp Files
%Requires Generic Bone inp Files
%Input frame
%Have Run TranformationMatrices.m
%Ensure if elseif statement for bone and body is correct

frame = 1;
chframe = num2str(frame);
ct_path = [path patient '/LtModel/'];
%creates a new folder for the frame of interest to put input files states in this frame
mkdir([ct_path name  'Frame' chframe]);

%load transformed node coordinates
if exist([trans_path name 'frame' num2str(frame) 'coord.mat'],'file')==2
    load([trans_path name 'frame' num2str(frame) 'coord.mat'])
load([trans_path name 'frame' num2str(frame) 'GScoord.mat']);
load([trans_path name 'frame' num2str(frame) 'Refcoord.mat']);

else
    
    Nodescord=struct;
    GSnodes=struct;
    REFnodes=struct;
    for h = 1:10
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
        elseif h == 5
            body = 'tibia';
            bone = 'fibula';
        elseif h == 6
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
        else
            break
        end
        
        switch body %select correct transformation matix
            case 'pelvis'; T = TM.T_pel(:,:,frame);
            case 'femur'; T = TM.T_fem(:,:,frame);
            case 'tibia'; T = TM.T_tib(:,:,frame);
            case 'patella'; T = TM.T_pat(:,:,frame);
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
           strncmp(bone, 'tibia',5)==1
        [nodes, gsnodes, refnodes]=bonereadjw(patient,trial,bone,path,T,frame);
        Nodescord.(bone)=nodes;
      GSnodes.(bone)= gsnodes;
       REFnodes.(bone)= refnodes;
        else
         nodes=bonereadjw(patient,trial,bone,path,T,frame);
        Nodescord.(bone)=nodes;
        end
        save([trans_path name 'frame' num2str(frame) 'coord.mat'],'Nodescord');
        save([trans_path name 'frame' num2str(frame) 'GScoord.mat'],'GSnodes');
                save([trans_path name 'frame' num2str(frame) 'Refcoord.mat'],'REFnodes');

    end
end
%%
disp('test4')
%create muscle connector elements for frame/s
frame=1;
muscleframecoordinates(path,patient,trial,frame);

%%
disp('test5')
%BONE PROPERTIES FILE -%dont need to modify since RBRN nodes already defined?
%print RBRN coord and beam to each other

fid=fopen([path patient '/LtModel/Props_Bone.inp'],'r'); %need generic file?
Z = textscan(fid, '%s', 'delimiter', '\n'); % import full file
z = [Z{1}]; % separate different lines into cells
fclose(fid);
format = '%-1.0f, %-1.6f, %-1.6f, %-1.6f';


%Bones ={  'pelvis',   'fibula'   'foot' 'bofoot' 'femoral_component' 'tibial_tray' }; %'patella' 
Bones ={  'pelvis',   'femur'  'tibia' 'patella'...
     'fibula'  'talus' 'foot' 'bofoot' 'femoral_component' 'tibial_tray' }; % 
for i = 1:size(Bones,2) %8 Bone Origins to be Defined
    bone=Bones{i};
    C=(Nodescord.(bone)(1,:)); %first node is RBRN
    %already labelled bones
%     switch i 
%         case 1; n= 800;
%         case 2; n=100;
%         case 3;  n=200;
%         case 4; n=300;  
%         case 5;  n=700;
%         case 6; n=900;
%         case 7;   n=1300;
%         case 8;  n=1400;
%         case 9;   n=400;
%         case 10;   n=500; %Select Bone Coordinates and Transformation Matrix
%     end
    %edit = [n C];
    Edit = num2str(C,format);
    z(i+1) = {Edit};
end
    fid = fopen([path patient '/LtModel/Props_Bone_' patient trial '.inp'],'w+');
fprintf(fid,'%s \n',z{:});
fclose(fid);
        %%
        %RBRN =rigid body ref node =num(i)
        %Beam to define a rigid beam connection to constrain the displacement and rotation of each slave node to the displacement and rotation of the control point.
        %         beam, 200, 700 Tib/fib
        % beam, 500, 1000
        
        % beam, 900, 1300 Foot/Tal
        % beam, 1300, 1400 foot/bofoot
        % beam, 1200, 1500
        % beam, 1500, 1600
        
        % Nodescord.(Bones(i));
        %c=node no?
        %Foot talus boffoot and fibula diffenrent mass element bone props
        %only need to do for Fib, pel,foot, bofoot -bones not in knee joint
%         case 1; T = TM.T_pel; num(i) = 800; 
%             %         case 2; T = T_fem; num(i) = 100; nnodes=92980;
%             %         case 3;  T = T_pat; num(i) = 300; nnodes=4617;
%             %         case 4; T = T_tib; num(i) = 200; nnodes=70322;
%             %        case 6;  T = T_tal; num(i) = 900; nnodes=2394;
%         case 2; T = TM.T_tib; num(i) = 700;   %fibular
%             
%         case 3;   T = TM.T_cal; num(i) = 1300;  %foot
%         case 4;  T = TM.T_toe; num(i) = 1400; %bofoot
%         case 5;   T=TM.T_FemComp; num(i) = 400; 
%         case 6;   T = TM.T_TTray; num(i) = 500; 
%             
%     end
%  
%     bone=char(Bones(i));
%        nnodes=size(Nodescord.(bone),1);
%     c=1; %Set RBRN as first node
%     file=([ct_path patient trial 'Frame' chframe '/' bone '_frame' chframe '.inp']); %inp file to write to
%     C = dlmread(file,',',[7 1 (nnodes+6) 3]); %node coordinates
%     %set RBRF to node no. num (i)
%     refs(i,1:4) = [1 C(c,:)];      refs(i,:) = T(:,:,1)*refs(i,:)';
%     edit = [num(i) refs(i,2:4)];
%     Edit = num2str(edit,format);
%     z(i+1) = {Edit};
% end
% 
% fid = fopen([path patient '/LtModel/Props_Bone_' patient trial '.inp'],'w+');
% fprintf(fid,'%s \n',z{:});
% fclose(fid);

%%
disp('test6')
% GS NODE AND KINEMATIC FILES

%GS Nodes and Kinematic Axes: TF, PF, Hip, and Ankle Joints of both legs
%Creates GS Nodes and GS Kinematic Axes inp files
%Have Run Transformation Maticies and Bone Coordinates
%Ensure correct selection of GS Nodes

%GS Nodes
clear refs
Bones ={    'femur' 'tibia' 'patella' 'pelvis',  'talus'  }; %
%define nodes as med/lat/inf/sup =num(i)+(1-4) node 00 is RBRN
clear refs
for i=[1:5]
    bone=Bones{i};
    for j=1:4
        refs((4*i+j-4),[1:4])=[1 GSnodes.(bone)(j,2:4)]; %gets untransformed GS node coordinates
       n((4*i+j-4),[1:4])=[1 REFnodes.(bone)(j,2:4)]; %already transformed
    end
end

%Reference Frame Nodes

%%
% n(1,1:4) = [1 Nodescord.('femur')(3,:)]; %Right Femur Origin
% n(2,1:4) = [1 (Nodescord.('femur')(3,:) + (TM.T_fem(:,:,1)*[1 1 0 0]')(2:4)];
% n(3,1:4) = [1 (Nodescord.('femur')(3,:) + [0 1 0])];
% n(4,1:4) = [1 (Nodescord.('femur')(3,:) + [0 0 1])];
% n(5,1:4) = [1 Nodescord.('tibia')(4,:)]; %Right Tibia Origin/sup
% n(6,1:4) = [1 (Nodescord.('tibia')(4,:) +[1 0 0])];
% n(7,1:4) = [1 (Nodescord.('tibia')(4,:) +[0 1 0])];
% n(8,1:4) = [1 (Nodescord.('tibia')(4,:) + [0 0 1])];
% n(9,1:4) = [1 Nodescord.('patella')(3,:)]; %Right Patella Origin/inf
% n(10,1:4) = [1 (Nodescord.('patella')(3,:) +[1 0 0])];
% n(11,1:4) = [1 (Nodescord.('patella')(3,:) +[0 1 0])];
% n(12,1:4) = [1 (Nodescord.('patella')(3,:) + [0 0 1])];
% n(13,1:4) = [1 Nodescord.('pelvis')(4,:)]; %Sacrum Origin/sup
% n(14,1:4) = [1 (Nodescord.('pelvis')(4,:) + [1 0 0])];
% n(15,1:4) = [1 (Nodescord.('pelvis')(4,:) + [0 1 0])];
% n(16,1:4) = [1 (Nodescord.('pelvis')(4,:) + [0 0 1])];
% n(17,1:4) = [1 Nodescord.('talus')(4,:)]; %Right Talus Origin/sup
% n(18,1:4) = [1 (Nodescord.('talus')(4,:) + [1 0 0])];
% n(19,1:4) = [1 (Nodescord.('talus')(4,:) + [0 1 0])];
% n(20,1:4) = [1 (Nodescord.('talus')(4,:) + [0 0 1])];
%%

%%
%ALREADY TRANSFORMED
% %Transformed Values and Node Numbers for origin nodes N=node no.s
% n(1,:)  = TM.T_fem(:,:,1)*n(1,:)'; N(1) = 1001;
% n(2,:)  = TM.T_fem(:,:,1)*n(2,:)'; N(2) = 1002;
% n(3,:)  = TM.T_fem(:,:,1)*n(3,:)'; N(3) = 1003;
% n(4,:)  = TM.T_fem(:,:,1)*n(4,:)'; N(4) = 1004;
% n(5,:)  = TM.T_tib(:,:,1)*n(5,:)'; N(5) = 2001;
% n(6,:)  = TM.T_tib(:,:,1)*n(6,:)'; N(6) = 2002;
% n(7,:)  = TM.T_tib(:,:,1)*n(7,:)'; N(7) = 2003;
% n(8,:)  = TM.T_tib(:,:,1)*n(8,:)'; N(8) = 2004;
% n(9,:)  = TM.T_pat(:,:,1)*n(9,:)'; N(9) = 3001;
% n(10,:)  = TM.T_pat(:,:,1)*n(10,:)'; N(10) = 3002;
% n(11,:)  = TM.T_pat(:,:,1)*n(11,:)'; N(11) = 3003;
% n(12,:)  = TM.T_pat(:,:,1)*n(12,:)'; N(12) = 3004;
% n(13,:)  = TM.T_pel(:,:,1)*n(13,:)'; N(13) = 8001; %changed from 1700
% n(14,:)  = TM.T_pel(:,:,1)*n(14,:)'; N(14) = 8002;
% n(15,:)  = TM.T_pel(:,:,1)*n(15,:)'; N(15) = 8003;
% n(16,:)  = TM.T_pel(:,:,1)*n(16,:)'; N(16) = 8004;
% n(17,:) = TM.T_tal(:,:,1)*n(17,:)'; N(17) = 9001;
% n(18,:) =TM.T_tal(:,:,1)*n(18,:)'; N(18) = 9002;
% n(19,:) = TM.T_tal(:,:,1)*n(19,:)'; N(19) = 9003;
% n(20,:) = TM.T_tal(:,:,1)*n(20,:)'; N(20) = 9004;
%%
%AP Nodes for all joints
%index for GS reference
for i = 1:4
    if i == 1 %Right TF Joint %tibia and femur
        j = 8; %Bone A Superior Node BONE A=tibia
        k = 7; %Bone A Inferior Node
        l = 2; %Bone B Lateral Node BONE B =femur
        m = 1; %Bone B Medial Node
        p = 8; %Bone A Ori Node
        o = 3; %Bone B Ori Node
    elseif i == 2 %Right PF Joint  %patella and femur 
        j = 12; %Bone A Superior Node A=patella B=femur
        k = 11; %Bone A Inferior Node
        l = 2; %Bone B Lateral Node
        m = 1; %Bone B Medial Node
        p = 11; %Bone A Ori Node
        o = 3; %Bone B Ori Node
        
    elseif i == 3 % Right Hip Joint %pelvis and femur
        j = 16; %Bone A Superior Node BONE A=HIP BONE B =femur
        k = 15; %Bone A Inferior Node
        l = 2; %Bone B Lateral Node
        m = 1; %Bone B Medial Node
        p = 14; %Bone A Ori Node
        o = 4; %Bone B Ori Node
        
    elseif i == 4 %Right Ankle Joint %talus and tibia
        j = 8; %Bone A Superior Node BONE A = tibia
        k = 7; %Bone A Inferior Node BONE B =talus
        l = 18; %Bone B Lateral Node
        m = 17; %Bone B Medial Node
        p = 7; %Bone A Ori Node
        o = 20; %Bone B Ori Node
        
    end
    
    % Joint descriptive axes
    % GS Joint: IE along Bone A axis, ML along Bone B axis, AP floating
    % (orthogonal to both) (Following Bull et al., KSSTA, 2002)
    X_ie = refs(j,2:4) - refs(k,2:4); X_ie = X_ie./sqrt(dot(X_ie,X_ie)); %Superior-Inferior Nodes Bone A
    X_ml = refs(l,2:4) - refs(m,2:4); X_ml = X_ml./sqrt(dot(X_ml,X_ml)); %Lateral-Medial Nodes Bone B
    X_ap = cross(X_ie,X_ml); X_ap = X_ap./sqrt(dot(X_ap,X_ap)); %Computationally determined AP axis
    
    % Axial plane
    plane_p = refs(o,2:4); % point on the plane : Bone B Oigin Node(Superior/Inferior whichever is closer to joint)
    plane_n = cross(X_ml,X_ap); % normal to the plane
    plane_n = plane_n./sqrt(dot(plane_n,plane_n));
    d = dot(plane_p,plane_n)*(-1);
    
    line_p = refs(p,2:4); % point on the line : Bone A Oigin Node(Superior/Inferior whichever is closer to joint)
    line_n = X_ie; % vector parallel to line
    % line parameter
    t = (-1)*(dot(plane_n,line_p) + d)/(dot(plane_n,line_n));
    x = line_p(1) + line_n(1)*t;
    y = line_p(2) + line_n(2)*t;
    z = line_p(3) + line_n(3)*t;
    inter_X1 = [x y z]; % AP NODE X50 : IE Axis with Bone A Node
    
    % Sagittal plane
    plane_p = refs(p,2:4); % point on the plane : Bone A Oigin Node
    plane_n = cross(X_ie,X_ap); % normal to the plane
    plane_n = plane_n./sqrt(dot(plane_n,plane_n));
    d = dot(plane_p,plane_n)*(-1);
    
    line_p = refs(o,2:4); % point on the line : Bone B Oigin Node
    line_n = X_ml; % vector parallel to line
    % line parameter
    t = (-1)*(dot(plane_n,line_p) + d)/(dot(plane_n,line_n));
    
    
    x = line_p(1) + line_n(1)*t;
    y = line_p(2) + line_n(2)*t;
    z = line_p(3) + line_n(3)*t;
    
    inter_X2 = [x y z]; % AP NODE X51 : ML Axis with Bone B Node
    
    if i == 1
        intertf1 = inter_X1;
        intertf2 = inter_X2;
    elseif i == 2
        interpf1 = inter_X1;
        interpf2 = inter_X2;
%     
    elseif i == 3
        interhip1 = inter_X1;
        interhip2 = inter_X2;
        
    elseif i == 4
        interank1 = inter_X1;
        interank2 = inter_X2;
        
    end
end


%%
disp('test7')
%Creates GS Nodes file -need to transform nodes? 

fout1=fopen([path patient '\LtModel\GS_Nodes_' name '.inp'],'w+');
fprintf(fout1,'*NODE, nset=gs_nodes');
for i = 1:20 %define GSnodes and coordinates =node no.s 1:20
    fprintf(fout1, '\n%1.0f, %1.5f, %1.5f, %1.5f', [i refs(i,2:4)]);
end

%node numbers for reference nodes. 
for i = 1:size(Bones,2) 
    switch i 
        case 1; I= 100;
        case 2; I=200;
        case 3;  I=300;
        case 4; I=800;  
        case 5;  I=900;
    end
    for j=1:4
        N(4*i+j-4)=[I+j]; 
        %
    end
end
fprintf(fout1,'\n**\n*mpc');
for i = 1:4
    fprintf(fout1,'\nbeam, '); %beams GS femur nodes nodes (1-4) to femur origin
    fprintf(fout1, '%1.0f, %1.0f', [i N(1)]);
end
for i = 5:8 %tibia -repeat beam for other bones
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [i N(5)]);
end
for i = 9:12 %patella
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [i N(9)]);
end
for i = 13:16 %pelvis
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [i N(13)]);
end
for i = 17:20 %talus
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [i N(17)]);
end

fprintf(fout1,'\n**');
fclose(fout1);


%%
disp('test8')
%Creates GS Kin Axes file

axes_list = {'R_TF_EXT_IE', 'R_TF_EXT_ML', 'R_TF_EXT_AP',  ...
    'R_PF_EXT_IE', 'R_PF_EXT_ML', 'R_PF_EXT_AP', ...
       'R_Hip_IE', 'R_Hip_ML', 'R_Hip_AP',  ...
    'R_Ank_IE', 'R_Ank_ML', 'R_Ank_AP'};
joint_list = {'R_TF_EXT', 'R_PF_EXT', 'R_Hip',  'R_Ank'};


fout1=fopen([path patient '\LtModel\GS_Kin_Axes_' name '.inp'],'w+');
fprintf(fout1,'**\n*NODE, nset=ref_frames\n');
for i = 1:20 %defines origin nodes
    fprintf(fout1,'%1.0f, %1.5f, %1.5f, %1.5f\n',[N(i) n(i,2:4)]);
end
fprintf(fout1,'**\n*mpc');
for i = 2:4
    fprintf(fout1,'\nbeam, '); %beams orgin reference nodes to each other for each bone
    fprintf(fout1, '%1.0f, %1.0f', [N(i) N(1)]);
end
for i = 6:8
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) N(5)]);
end
for i = 10:12
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) N(9)]);
end
for i = 14:16
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) N(13)]);
end
for i = 18:20
    fprintf(fout1,'\nbeam, ');
    fprintf(fout1, '%1.0f, %1.0f', [N(i) N(17)]);
end


for i = 1:4
    fprintf(fout1,['\n**\n**\n**Intersection nodes to define floating AP axes - ' joint_list{i} '\n*NODE\n']);
    switch joint_list{i} %prints node no.s and coordinates
        case 'R_TF_EXT'; Node = [150 151]; Cor1 = intertf1; Cor2 = intertf2;
       case 'R_PF_EXT'; Node = [250 251]; Cor1 = interpf1; Cor2 = interpf2;
        
        case 'R_Hip'; Node = [350 351]; Cor1 = interhip1; Cor2 = interhip2;
        
        case 'R_Ank'; Node = [450 451]; Cor1 = interank1; Cor2 = interank2;
        
    end
    fprintf(fout1, '%1.0f, %1.5f, %1.5f, %1.5f', [Node(1) Cor1]);
    fprintf(fout1, '\n%1.0f, %1.5f, %1.5f, %1.5f', [Node(2) Cor2]);
end

for i = 1:4
    switch joint_list{i}
        case 'R_TF_EXT'; Node = [150 151]; Element = [150 151 152]; Ori = [2001 1001]; %origin nodes
        case 'R_PF_EXT'; Node = [250 251]; Element = [250 251 252]; Ori = [300 100];
        
        case 'R_Hip'; Node = [350 351]; Element = [350 351 352]; Ori = [8001 4]; %4=GS node femur hip/superior
        
        case 'R_Ank'; Node = [450 451]; Element = [450 451 452]; Ori = [7 9001]; %7=GS tibiba  inferior node
        
    end
    fprintf(fout1,['\n**\n**\n********** ' joint_list{i} ' Kinematic Axes\n**\n*ELEMENT,TYPE=CONN3D2,ELSET=axis_' joint_list{i} '_IE\n']);
    fprintf(fout1, '%1.0f, %1.0f, %1.0f', [Element(1) Ori(1) Node(1)]);
    fprintf(fout1,['\n*ELEMENT,TYPE=CONN3D2,ELSET=axis_' joint_list{i} '_ML\n']);
    fprintf(fout1, '%1.0f, %1.0f, %1.0f', [Element(3) Node(2) Ori(2)]);
    fprintf(fout1,['\n*ELEMENT,TYPE=CONN3D2,ELSET=axis_' joint_list{i} '_AP\n']);
    fprintf(fout1, '%1.0f, %1.0f, %1.0f', [Element(2) Node(1) Node(2)]);
    
    
    if i == 1  %TF Joints
        for j = 1:3
            k = 3*(i-1)+j;
            fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
            if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(1) Ori(1)+1 Ori(1)]);
            elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(2) Ori(2)+2 Node(2)]);
            elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(2) Ori(2) Node(1)]);
            end
        end
    elseif i == 2 
        for j = 1:3
            k = 3*(i-1)+j;
            fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
            if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(1) Ori(1)+1 Node(1)]);
            elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(2) Ori(2)+2 Ori(2)]);
            elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(1) Ori(2) Node(2)]);
            end
        end
    elseif i == 3  %Hip Joints
        for j = 1:3
            k = 3*(i-1)+j;
            fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
            if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(1) Ori(1)+1 Node(1)]);
            elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(2) Ori(2)+2 Node(2)]);
            elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(2) Ori(2) Node(1)]);
            end
        end
    elseif i == 4  %Ankle Joints
        for j = 1:3
            k = 3*(i-1)+j;
            fprintf(fout1,['\n**\n**\n*CONNECTOR SECTION, ELSET=axis_' axes_list{k} '\nCYLINDRICAL,\naxis_' axes_list{k} '_ORI,\n*ORIENTATION, NAME=axis_' axes_list{k} '_ORI, DEFINITION=NODES']);
            if j == 1; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(1) Ori(1)+1 Node(1)]);
            elseif j == 2; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Ori(2) Ori(2)+2 Node(2)]);
            elseif j ==3; fprintf(fout1, '\n%1.0f, %1.0f, %1.0f', [Node(1) Ori(2) Node(2)]);
            end
        end
    end
   
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

disp('test10.5') 
% 
 Musc_amp_list = {'addbrev' 'addlong' 'addmagProx' 'addmagMid' 'addmagDist' 'addmagIsch'...
    'bflh' 'bfsh' 'edl' 'ehl' 'fdl' 'fhl' 'gaslat' 'gasmed' 'gem'...
    'glmax1' 'glmax2' 'glmax3' 'glmed1' 'glmed2' 'glmed3' 'glmin1' 'glmin2' 'glmin3' ...
    'grac' 'iliacus' 'pect' 'perbrev' 'perlong' 'pertert' 'piri' 'psoas' ...
    'quadfem' 'recfem' 'sart' 'semimem' 'semiten' 'soleus' 'tfl' 'tibant'...
    'tibpost' 'vasint' 'vaslat' 'vasmed'};

%creates muscle amplitude files for trial
muscfile=[amp_path patient  trial '-MUSC-BL.inp'];
if exist(muscfile,'file') ~= 2
    disp('creating amplitute card')
    MUSC=MUSC_AMP_CARDjw(path,patient,trial);
end
%apply amplitudes to connectors
fout1=fopen([path patient '\MUSC\apply_amp_musc' name '_R.inp'],'w+');
for var=1:size(Musc_amp_list,2) 
 NAME = Musc_amp_list{var};
    fprintf(fout1,['*CONNECTOR LOAD, OP=NEW, AMPLITUDE=' NAME '_amp\n']);
    fprintf(fout1,['MUSC_' NAME ', 1, 0\n']); %need O baseload??
    
      end
fprintf(fout1,'**\n');


%%
disp('test17') %like test 11
% AMPLITUDE CARDS
%Initial Flexion -get into starting position -already is so dont need
%%
% %Initial Loads - Hypermesh Knee
% %file = [path patient '\GEOM\BONE-TIBIA_R-' name '.inp'];
% coorTibia = Nodescord.('tibia'); %dlmread(file,',',[1 0 15272 3]);
% % file = [path patient '\MUSC\GS_KIN_AXIS_R-' name '.inp'];
% % coorAxis = dlmread(file,',',[3 0 18 3]);
% file = [path patient '\LtModel\GS_Nodes_' name '.inp'];
% coorGS = dlmread(file,',',[1 0 16 3]);
% 
% gsF = [1 coorGS(1,2:4); 1 coorGS(2,2:4); 1 coorGS(3,2:4); 1 coorGS(4,2:4)];
% gsT = [1 coorGS(5,2:4); 1 coorGS(6,2:4); 1 coorGS(7,2:4); 1 coorGS(8,2:4)];
% Ts = coorGS(8,2:4); %superior node of tibia
% 
% 
% %Bone A IE HM-Tibia; Bone B ML LT-Femur
% % Joint descriptive axes
% % GS Joint: IE along Bone A axis, ML along Bone B axis, AP floating
% % (orthogonal to both) (Following Bull et al., KSSTA, 2002)
% X_si = gsT(4,2:4) - gsT(3,2:4); X_si = X_si./sqrt(dot(X_si,X_si)); %Superior-Inferior Nodes Bone A
% X_ml = gsF(2,2:4) - gsF(1,2:4); X_ml = X_ml./sqrt(dot(X_ml,X_ml)); %Lateral-Medial Nodes Bone B
% X_ap = cross(X_si,X_ml); X_ap = X_ap./sqrt(dot(X_ap,X_ap)); %Computationally determined AP axis
% 
% % Axial plane
% plane_p = gsF(3,2:4); % point on the plane : Bone B Oigin Node(Superior/Inferior whichever is closer to joint)
% plane_n = cross(X_ml,X_ap); % normal to the plane
% plane_n = plane_n./sqrt(dot(plane_n,plane_n));
% d = dot(plane_p,plane_n)*(-1);
% 
% line_p = gsT(4,2:4); % point on the line : Bone A Oigin Node(Superior/Inferior whichever is closer to joint)
% line_n = X_si; % vector parallel to line
% % line parameter
% t = (-1)*(dot(plane_n,line_p) + d)/(dot(plane_n,line_n));
% x = line_p(1) + line_n(1)*t; y = line_p(2) + line_n(2)*t; z = line_p(3) + line_n(3)*t;
% inter_X1 = [x y z]; % AP NODE X50 : IE Axis with Bone A Node
% 
% % Sagittal plane
% plane_p = gsT(4,2:4); % point on the plane : Bone A Oigin Node
% plane_n = cross(X_si,X_ap); % normal to the plane
% plane_n = plane_n./sqrt(dot(plane_n,plane_n));
% d = dot(plane_p,plane_n)*(-1);
% 
% line_p = gsF(3,2:4); % point on the line : Bone B Oigin Node
% line_n = X_ml; % vector parallel to line
% % line parameter
% t = (-1)*(dot(plane_n,line_p) + d)/(dot(plane_n,line_n));
% 
% x = line_p(1) + line_n(1)*t;
% y = line_p(2) + line_n(2)*t;
% z = line_p(3) + line_n(3)*t;
% 
% inter_X2 = [x y z]; % AP NODE X51 : ML Axis with Bone B Node
% 
% hmlT = gsT(2,2:4) - gsT(1,2:4); hmlT = hmlT./sqrt(dot(hmlT,hmlT)); %temporary ML axis of the tibia
% hsi = gsT(4,2:4) - gsT(3,2:4);hsi = hsi./sqrt(dot(hsi,hsi))*(-1); %SI axis of the tibia, used as the SI axis of the joint
% hap = cross(hmlT,hsi); hap = hap./sqrt(dot(hap,hap)); %AP axis of the tibia
% hml = cross(hsi,hap); hml = hml./sqrt(dot(hml,hml)); %final ML axis of the tibia
% H1=zeros(4,4); H1(1,1)=1; H1(2:4,1)=gsT(4,2:4)'; H1(2:4,2)=hml'; H1(2:4,3)=hap'; H1(2:4,4)=hsi'; %Matrix of the tibia axes and point of origin(inferior point of the pelvis)
% 
% % Femoral Reference Frame
% fmlT = gsF(2,2:4) - gsF(1,2:4); fmlT = fmlT./sqrt(dot(fmlT,fmlT));
% fsi = gsF(4,2:4) - gsF(3,2:4); fsi = fsi./sqrt(dot(fsi,fsi))*(-1);
% fap = cross(fmlT,fsi); fap = fap./sqrt(dot(fap,fap));
% fml = cross(fsi,fap); fml = fml./sqrt(dot(fml,fml));
% F1=zeros(4,4); F1(1,1)=1; F1(2:4,1)=gsF(3,2:4)'; F1(2:4,2)=fml'; F1(2:4,3)=fap'; F1(2:4,4)=fsi'; %origin of the femur is the superior point
% 
% RHF=F1\H1; %Transformation from local femur into local pelvis
% 
% %TF kinematics: formulas based on Grood-Suntay paper
% IkinH(1,1)=1*atan(RHF(3,4)/RHF(4,4)); %FE angle
% if IkinH(1,1)<-.5 %To flip sign when trialon goes past 90 degrees
%     IkinH(1,1)=IkinH(1,1)+pi;
%     disp('FE flipped')
% end
% IkinH(1,2)=acos(RHF(2,4))-pi/2; %VV angle
% IkinH(1,2)=-IkinH(1,2);
% IkinH(1,3)=atan(RHF(2,3)/RHF(2,2)); %IE angle
% IkinH(1,3)=-IkinH(1,3);
% IkinH(1,4)=RHF(2,1); %ML translation
% IkinH(1,5)=RHF(3,1)*cos(IkinH(1,1))-RHF(4,1)*sin(IkinH(1,1)); %AP translation (text has TF(1,4))
% IkinH(1,6)=-(RHF(2,4)*RHF(2,1)+RHF(3,4)*RHF(3,1)+RHF(4,4)*RHF(4,1)); %SI translation

% %Initial Loads - LT knee
% %file = ['LtModel\' name 'Frame1\femur_r_frame1.inp'];
% file = [path patient '\LtModel\' name 'Frame1\femur_frame1.inp'];
% coorHMF = dlmread(file,',',[7 0 462 3]);
% file = [path patient '\LtModel\' name 'Frame1\tibia_frame1.inp'];
% coorTibia = dlmread(file,',',[7 0 151 3]);
% % file = ['MUSC\GS_KIN_AXIS-' name '.inp'];
% % coorAxis = dlmread(file,',',[3 0 3 3]);
% 
% gsF = [1 coorHMF(434,2:4); 1 coorHMF(411,2:4); 1 coorHMF(421,2:4); 1 coorHMF(2,2:4)];
% gsT = [1 coorTibia(17,2:4); 1 coorTibia(2,2:4); 1 coorTibia(145,2:4); 1 coorTibia(1,2:4)];
% %Femur m:4020, l:17836, i:400000, s:10074
% %Tibia m:45262, l:55368, i:52148, s:50535
% 
% hmlT = gsT(2,2:4) - gsT(1,2:4); hmlT = hmlT./sqrt(dot(hmlT,hmlT)); %temporary ML axis of the tibia
% hsi = gsT(4,2:4) - gsT(3,2:4);hsi = hsi./sqrt(dot(hsi,hsi))*(-1); %SI axis of the tibia, used as the SI axis of the joint
% hap = cross(hmlT,hsi); hap = hap./sqrt(dot(hap,hap)); %AP axis of the tibia
% hml = cross(hsi,hap); hml = hml./sqrt(dot(hml,hml)); %final ML axis of the tibia
% H1=zeros(4,4); H1(1,1)=1; H1(2:4,1)=gsT(4,2:4)'; H1(2:4,2)=hml'; H1(2:4,3)=hap'; H1(2:4,4)=hsi'; %Matrix of the tibia axes and point of origin(inferior point of the pelvis)
% 
% % Femoral Reference Frame
% fmlT = gsF(2,2:4) - gsF(1,2:4); fmlT = fmlT./sqrt(dot(fmlT,fmlT));
% fsi = gsF(4,2:4) - gsF(3,2:4); fsi = fsi./sqrt(dot(fsi,fsi))*(-1);
% fap = cross(fmlT,fsi); fap = fap./sqrt(dot(fap,fap));
% fml = cross(fsi,fap); fml = fml./sqrt(dot(fml,fml));
% F1=zeros(4,4); F1(1,1)=1; F1(2:4,1)=gsF(3,2:4)'; F1(2:4,2)=fml'; F1(2:4,3)=fap'; F1(2:4,4)=fsi'; %origin of the femur is the inferior point
% 
% RHF=inv(F1)*H1; %Transformation from local femur into local pelvis
% 
% %TF kinematics: formulas based on Grood-Suntay paper
% IkinL(1,1)=1*atan(RHF(3,4)/RHF(4,4)); %FE angle
% if IkinL(1,1)<-.5 %To flip sign when trialon goes past 90 degrees
%     IkinL(1,1)=IkinL(1,1)+pi;
%     disp('FE flipped')
% end
% IkinL(1,2)=acos(RHF(2,4))-pi/2; %VV angle
% IkinL(1,3)=atan(RHF(2,3)/RHF(2,2)); %IE angle
% IkinL(1,4)=RHF(2,1); %ML translation
% IkinL(1,5)=(RHF(3,1)*cos(IkinL(1,1)))-(RHF(4,1)*sin(IkinL(1,1))); %AP translation (text has TF(1,4))
% IkinL(1,6)=-(RHF(2,4)*RHF(2,1)+RHF(3,4)*RHF(3,1)+RHF(4,4)*RHF(4,1)); %SI translation

%%
%Write Initial Amp Card test 17 continued for both kinematics and muscles
fout1=fopen([path patient '\MUSC\amp_profiles_initial' name '_R.inp'],'w+');

t = [0.0 1]; %changes from 0.57

%%
% load_values = [-0.1667, -0.1111, -0.1, -0.0429, -0.12, -0.0571];
% LPFL = '-5.0';
% MPFL = '-16.0';

% 
% f1 = zeros(2,6); %f2 = zeros(2,5);
% f1(1,1:6) = IkinH(1:6);
% %f1(2,1:6) = IkinL(1:6);
% %f1(:,4) = [ML;ml]; f1(:,5) = [AP;ap]; f1(:,6) = [SI,si];
% %f2(2,1) = musc_R(1,1); f2(2,2) = musc_R(1,2); f2(2,3) = musc_R(1,3); f2(2,4) = musc_R(1,4); f2(2,5) = 1;
% 
% %kinematics
% amp_list = {'AMP_RTF_FE_LOAD'; 'AMP_RTF_VV_LOAD'; 'AMP_RTF_IE_LOAD'; 'AMP_RTF_ML_LOAD'; 'AMP_RTF_AP_LOAD'; 'AMP_RTF_SI_LOAD'};
% 
% for var = 1:6
%     NAME = amp_list{var};
%     fprintf(fout1,['*AMPLITUDE, VALUE=RELATIVE, DEFINITION=SMOOTH STEP, NAME=' NAME '\n']);
%     %fprintf(fout1,'%1.4f,%1.4f,%1.4f,%1.4f,%1.4f,%1.4f\n',[t(1) 0.0000 0.5*t(2) 0.1*(f1(2,var)-f1(1,var)) t(2) f1(2,var)-f1(1,var)]);
%     fprintf(fout1,'%1.4f,%1.4f,%1.4f,%1.4f,%1.4f,%1.4f\n',[t(1) 0.0000 0.5*t(2) 0.1*(f1(1,var)) t(2) f1(1,var)]);
% 
% end
% fprintf(fout1,'**\n');

%%
% %muscle forces -initial amplitude card- need to remove initial ramp up
% when creating card so that get initial value. 

for var=1:size(Musc_amp_list,2) 
 NAME = Musc_amp_list{var};
    fprintf(fout1,['*AMPLITUDE, VALUE=RELATIVE, DEFINITION=SMOOTH STEP, NAME=' NAME '_int_amp\n']);
        %fprintf(fout1,'%1.4f,%1.4f,%1.4f,%1.4f,%1.4f,%1.4f\n',[t(1) 0.0000 0.5*t(2) 0.1*(f2(2,var)) t(2) f2(2,var)]);

    fprintf(fout1,'%1.4f,%1.4f,%1.4f,%1.4f,%1.4f,%1.4f,%1.4f,%1.4f\n%1.4f,%1.4f,%1.4f,%1.4f\n',[t(1) 0.0*MUSC(1,var) 0.1*t(2) 0.0*MUSC(1,var) 0.2*t(2) 0.001*MUSC(1,var) 0.5*t(2) 0.01*MUSC(1,var) 0.75*t(2) 0.1*MUSC(1,var) t(2) MUSC(1,var)]);
end
fprintf(fout1,'**\n');

% load_list = {'RF_LOAD = '; 'VI_LOAD = '; 'VLL_LOAD = '; 'VLO_LOAD = '; 'VML_LOAD = '; 'VMO_LOAD = '; 'LPFL_LOAD = '; 'MPFL_LOAD = '};
% fprintf(fout1,'*PARAMETER\n');
% for i = 1:6
%     fprintf(fout1,[load_list{i}]);
%     fprintf(fout1,'%1.4f\n',load_values(i)*musc_multi);
% end
% fprintf(fout1,[load_list{7} LPFL '\n']);
% fprintf(fout1,[load_list{8} MPFL '\n']);

fclose(fout1);

fout1=fopen([path patient '\MUSC\apply_amp_musc_initial' name '_R.inp'],'w+');
for var=1:size(Musc_amp_list,2) 
 NAME = Musc_amp_list{var};
    fprintf(fout1,['*CONNECTOR LOAD, OP=NEW, AMPLITUDE=' NAME '_int_amp\n']);
    fprintf(fout1,['MUSC_' NAME ', 1, 0\n']); %need O baseload??
    
      end
fprintf(fout1,'**\n');

fclose(fout1);
%Modify HyperMesh GS Kin Axes -dont need?
%%
% %open generic file and create cell array of contents
% clear B
% fid = fopen([path patient '\MUSC\GS_KIN_AXIS_R-' name '.inp']);
% tline = fgetl(fid);
% i = 0;
% while ischar(tline)
%     i = i+1;
%     tline = fgetl(fid);
%     B{i} = tline;
% end
% fclose(fid);
% 
% %edit generic cell array and export to a created frame file
% formatSpec = '%3.0f, %5.5f, %5.5f, %5.5f';
% %need to add a 0 since increased node numbers. 
% B{13} = sprintf(formatSpec,[10400004 Ts]); %tibia superior node
% B{17} = sprintf(formatSpec,[10400008 inter_X1]); % AP NODE X50 : IE Axis with Bone A Node
% B{18} = sprintf(formatSpec,[10400009 inter_X2]);  % AP NODE X51 : ML Axis with Bone B Node
% 
% 
% fid = fopen([path patient '\MUSC\GS_KIN_AXIS_R-' name 'Modify' version '.inp'], 'w+');
% for i = 1:numel(B)
%     if B{i+1} == -1
%         fprintf(fid,'%s', B{i});
%         break
%     else
%         fprintf(fid,'%s\n', B{i});
%     end
% end
% fclose(fid);

%%
disp('test13')
%Transform coordinates to global coordinate system of each frame; combine into one matrix
%to create kinematic amplitude card for gait cycle
gs = zeros(size(TM.T_tib,3),61); %=no. GS nodes*3 -1
gs(:,1) = 1;

for i=1:5
    bone=Bones{i};
    switch i %Select Bone Coordinates and Transformation Matrix
        case 1; T = TM.T_fem;
        case 2; T = TM.T_tib;
        case 3;   T = TM.T_pat;
        case 4;  T = TM.T_pel;
        case 5;   T=TM.T_tal;
    end
    
    %                 for j=1:4
    %                    GSc.(bone)(j,:)=[1 Nodescord.(bone)(j+1,:)]
    %                 end
    
    for k = 1:size(T,3) %num frames
        for j=1:4
            %for each GS node, and for each frame, multiply coordinates by transformation matrix
            GSTrans.(bone)(j,:)=T(:,:,k)*([1 GSnodes.(bone)(j,2:4)])';
        end
        gs(k,12*i-10:12*i+1)=[GSTrans.(bone)(1,2:4)  GSTrans.(bone)(2,2:4) GSTrans.(bone)(3,2:4) GSTrans.(bone)(4,2:4)];
    end
end

%AP Nodes for all joints
for i = 1:size(TM.T_tib,3)
    for a = 1:4
        %a gives joint no.
        %j-o gives GS ref node no.
        if a == 1                           %Right TF Joint
            j = 8;  %Bone A Superior Node
            k = 7;  %Bone A Inferior Node
            l = 2;  %Bone B Lateral Node
            m = 1;  %Bone B Medial Node
            p = 8;  %Bone A Ori Node
            o = 3;  %Bone B Ori Node
        elseif a == 2                       %Right PF Joint
            j = 12; %Bone A Superior Node
            k = 11; %Bone A Inferior Node
            l = 2;  %Bone B Lateral Node
            m = 1;  %Bone B Medial Node
            p = 11; %Bone A Ori Node
            o = 3;  %Bone B Ori Node
        
        elseif a == 3                       % Right Hip Joint
            j = 16; %Bone A Superior Node
            k = 15; %Bone A Inferior Node
            l = 2;  %Bone B Lateral Node
            m = 1;  %Bone B Medial Node
            p = 16; %Bone A Ori Node
            o = 4;  %Bone B Ori Node
      
        elseif a == 4                       %Right Ankle Joint
            j = 8;  %Bone A Superior Node
            k = 7;  %Bone A Inferior Node
            l = 18; %Bone B Lateral Node
            m = 17; %Bone B Medial Node
            p = 7;  %Bone A Ori Node
            o = 20; %Bone B Ori Node
       
        end
        
        % Joint descriptive axes
        % GS Joint: IE along Bone A axis, ML along Bone B axis, AP floating
        % (orthogonal to both) (Following Bull et al., KSSTA, 2002)
        X_si = gs(i,(3*j)-1:(3*j)+1) - gs(i,(3*k)-1:(3*k)+1); X_si = X_si./sqrt(dot(X_si,X_si)); %Superior-Inferior Nodes Bone A
        X_ml = gs(i,(3*l)-1:(3*l)+1) - gs(i,(3*m)-1:(3*m)+1); X_ml = X_ml./sqrt(dot(X_ml,X_ml)); %Lateral-Medial Nodes Bone B
        X_ap = cross(X_si,X_ml); X_ap = X_ap./sqrt(dot(X_ap,X_ap)); %Computationally determined AP axis
        
        % Axial plane
        plane_p = gs(i,(3*o)-1:(3*o)+1); % point on the plane : Bone B Oigin Node(Superior/Inferior whichever is closer to joint)
        plane_n = cross(X_ml,X_ap); % normal to the plane
        plane_n = plane_n./sqrt(dot(plane_n,plane_n));
        d = dot(plane_p,plane_n)*(-1);
        
        line_p = gs(i,(3*p)-1:(3*p)+1); % point on the line : Bone A Oigin Node(Superior/Inferior whichever is closer to joint)
        line_n = X_si; % vector parallel to line
        % line parameter
        t = (-1)*(dot(plane_n,line_p) + d)/(dot(plane_n,line_n));
        x = line_p(1) + line_n(1)*t; y = line_p(2) + line_n(2)*t; z = line_p(3) + line_n(3)*t;
        inter_X1 = [x y z]; % AP NODE X50 : IE Axis with Bone A Node
        
        % Sagittal plane
        plane_p = gs(i,(3*p)-1:(3*p)+1); % point on the plane : Bone A Oigin Node
        plane_n = cross(X_si,X_ap); % normal to the plane
        plane_n = plane_n./sqrt(dot(plane_n,plane_n));
        d = dot(plane_p,plane_n)*(-1);
        
        line_p = gs(i,(3*o)-1:(3*o)+1); % point on the line : Bone B Oigin Node
        line_n = X_ml; % vector parallel to line
        % line parameter
        t = (-1)*(dot(plane_n,line_p) + d)/(dot(plane_n,line_n));
        
        x = line_p(1) + line_n(1)*t;
        y = line_p(2) + line_n(2)*t;
        z = line_p(3) + line_n(3)*t;
        
        inter_X2 = [x y z]; % AP NODE X51 : ML Axis with Bone B Node
        
        if a == 1
            inter_rtf1(i,:) = inter_X1;
            inter_rtf2(i,:) = inter_X2;
        elseif a == 2
            inter_rpf1(i,:) = inter_X1;
            inter_rpf2(i,:) = inter_X2;
       
        elseif a == 3
            inter_rhip1(i,:) = inter_X1;
            inter_rhip2(i,:) = inter_X2;
        
        elseif a == 4
            inter_rank1(i,:) = inter_X1;
            inter_rank2(i,:) = inter_X2;
        
        end
    end
end       

kin=zeros(size(gs,1),24);
for i = 1:size(gs,1)
    gsRF = [1 gs(i,2:4); 1 gs(i,5:7); 1 gs(i,8:10); 1 gs(i,11:13)]  ; %RF m,l,i,s
    gsRT = [1 gs(i,14:16); 1 gs(i,17:19); 1 gs(i,20:22); 1 gs(i,23:25)]  ; %RT m,l,i,s
    gsRP = [1 gs(i,26:28); 1 gs(i,29:31); 1 gs(i,32:34); 1 gs(i,35:37)]  ; %RP m,l,i,s
    gsRH = [1 gs(i,41:43); 1 gs(i,38:40); 1 gs(i,44:46); 1 gs(i,47:49)]  ; %RHip m,l,i,s
    gsRA = [1 gs(i,50:52); 1 gs(i,53:55); 1 gs(i,56:58); 1 gs(i,59:61)]  ; %RAnk m,l,i,s
   
    % GS kinematics from GS probed points
    % right femoral reference frame
    rfmlT = gsRF(2,2:4) - gsRF(1,2:4); %temporary ML axis: lateral - medial
    rfmlT = rfmlT./sqrt(dot(rfmlT,rfmlT));
    
    rfsi = gsRF(4,2:4) - gsRF(3,2:4); %IE axis: superior - inferior
    rfsi = rfsi./sqrt(dot(rfsi,rfsi));
    
    rfap = cross(rfmlT,rfsi); %AP axis
    rfap = rfap./sqrt(dot(rfap,rfap));
    
    rfml = cross(rfsi,rfap); %ML axis
    rfml = rfml./sqrt(dot(rfml,rfml));
    RF1=zeros(4,4); RF1(1,1)=1; RF1(2:4,1)=gsRF(3,2:4)'; RF1(2:4,2)=rfml'; RF1(2:4,3)=rfap'; RF1(2:4,4)=rfsi';
    RF2=RF1; RF2(2:4,1)=gsRF(4,2:4)'; %Switched origin to superior node for hip joint
    
    
    % right tibial reference frame
    rtmlT = gsRT(2,2:4) - gsRT(1,2:4); %temporary ML axis: lateral - medial
    rtmlT = rtmlT./sqrt(dot(rtmlT,rtmlT));
    
    rtsi = gsRT(4,2:4) - gsRT(3,2:4); %IE axis: superior - inferior
    rtsi = rtsi./sqrt(dot(rtsi,rtsi));
    
    rtap = cross(rtmlT,rtsi); %AP axis
    rtap = rtap./sqrt(dot(rtap,rtap));
    
    rtml = cross(rtsi,rtap); %ML axis
    rtml = rtml./sqrt(dot(rtml,rtml));
    RT1=zeros(4,4); RT1(1,1)=1; RT1(2:4,1)=gsRT(4,2:4)'; RT1(2:4,2)=rtml'; RT1(2:4,3)=rtap'; RT1(2:4,4)=rtsi';
    RT2=RT1; RT2(2:4,1)=gsRT(3,2:4)'; %Switched origin to inferior node for ankle joint
    
    
    % right patella reference frame
    rpmlT = gsRP(2,2:4) - gsRP(1,2:4); %temporary ML axis: lateral - medial
    rpmlT = rpmlT./sqrt(dot(rpmlT,rpmlT));
    
    rpsi = gsRP(4,2:4) - gsRP(3,2:4); %IE axis: superior - inferior
    rpsi = rpsi./sqrt(dot(rpsi,rpsi));
    
    rpap = cross(rpsi,rpmlT); %AP axis
    rpap = rpap./sqrt(dot(rpap,rpap));
    
    rpml = cross(rpap,rpsi); %ML axis
    rpml = rpml./sqrt(dot(rpml,rpml));
    RP1=zeros(4,4); RP1(1,1)=1; RP1(2:4,1)=gsRP(3,2:4)'; RP1(2:4,2)=rpml'; RP1(2:4,3)=rpap'; RP1(2:4,4)=rpsi';
    
    
    % right pelvis/sacrum reference frame
    rhmlT = gsRH(2,2:4) - gsRH(1,2:4); %temporary ML axis: lateral - medial
    rhmlT = rhmlT./sqrt(dot(rhmlT,rhmlT));
    
    rhsi = gsRH(4,2:4) - gsRH(3,2:4); %IE axis: superior - inferior
    rhsi = rhsi./sqrt(dot(rhsi,rhsi));
    
    rhap = cross(rhsi,rhmlT); %AP axis
    rhap = rhap./sqrt(dot(rhap,rhap));
    
    rhml = cross(rhap,rhsi); %ML axis
    rhml = rhml./sqrt(dot(rhml,rhml));
    RH1=zeros(4,4); RH1(1,1)=1; RH1(2:4,1)=gsRH(4,2:4)'; RH1(2:4,2)=rhml'; RH1(2:4,3)=rhap'; RH1(2:4,4)=rhsi';
    
    
    % right talus reference frame
    ramlT = gsRA(2,2:4) - gsRA(1,2:4); %temporary ML axis: lateral - medial
    ramlT = ramlT./sqrt(dot(ramlT,ramlT));
    
    rasi = gsRA(4,2:4) - gsRA(3,2:4); %IE axis: superior - inferior
    rasi = rasi./sqrt(dot(rasi,rasi));
    
    raap = cross(rasi,ramlT); %AP axis
    raap = raap./sqrt(dot(raap,raap));
    
    raml = cross(raap,rasi); %ML axis
    raml = raml./sqrt(dot(raml,raml));
    RA1=zeros(4,4); RA1(1,1)=1; RA1(2:4,1)=gsRA(4,2:4)'; RA1(2:4,2)=raml'; RA1(2:4,3)=raap'; RA1(2:4,4)=rasi';
    
    
    
    RFT=inv(RT1)*RF1; %Transformation from local right tibia into local right femur
    RFP=inv(RP1)*RF1; %Transformation from local right patella into local right femur
    RHF=inv(RF2)*RH1;  %Transformation from local right femur into local right pelvis/sacrum
    RTA=inv(RA1)*RT2; %Transformation from local right talus into local right tibia
    

    %Right TF kinematics
    kin(i,1)=1*atan(RFT(3,4)/RFT(4,4)); %FE angle
    if kin(i,1)>0 %To flip sign when trialon goes past 90 degrees
        kin(i,1)=kin(i,1)-pi;
    end
    kin(i,2)=acos(RFT(2,4))-pi/2; %VV angle
    kin(i,3)=atan(RFT(2,3)/RFT(2,2)); %IE angle
    kin(i,4)=RFT(2,1); %ML translation
    kin(i,5)=RFT(3,1)*cos(kin(i,1))-RFT(4,1)*sin(kin(i,1)); %AP translation (text has TF(1,4))
    kin(i,6)=(RFT(2,4)*RFT(2,1)+RFT(3,4)*RFT(3,1)+RFT(4,4)*RFT(4,1)); %SI translation
    
    p1 = gs(i,23:25); p2 = inter_rtf1(i,:);
    dist = sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2 + (p1(3) - p2(3))^2);
    deltaML = sin(kin(i,2))*dist;
    kin(i,4) = kin(i,4) + deltaML;
    
    
    %Right PF kinematics
    kin(i,7)=1*atan(RFP(3,4)/RFP(4,4)); %patFE angle
    %        if kin(trial,7)<-.5 %To flip sign when trialon goes past 90 degrees
    %            kin(trial,7)=kin(trial,7)+pi;
    %        end
    kin(i,8)=acos(RFP(2,4))-pi/2; %patVV angle
    kin(i,9)=atan(RFP(2,3)/RFP(2,2)); %patIE angle
    kin(i,10)=RFP(2,1); %patML translation
    kin(i,11)=RFP(3,1)*cos(kin(i,7))-RFP(4,1)*sin(kin(i,7)); %patAP translation (text has TF(1,4))
    kin(i,12)=(RFP(2,4)*RFP(2,1)+RFP(3,4)*RFP(3,1)+RFP(4,4)*RFP(4,1)); %patSI translation
    
    p1 = gs(i,32:34); p2 = inter_rpf1(i,:);
    dist = sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2 + (p1(3) - p2(3))^2);
    deltaML = sin(kin(i,8))*dist;
    kin(i,10) = kin(i,10) + deltaML;
    
    
    %Right Hip kinematics
    kin(i,13)=1*atan(RHF(3,4)/RHF(4,4)); %FE angle
    if kin(i,13)>0 %To flip sign when trialon goes past 90 degrees
        kin(i,13)=kin(i,13)-pi;
    end
    kin(i,14)=acos(RHF(2,4))-pi/2; %VV angle
    kin(i,15)=atan(RHF(2,3)/RHF(2,2)); %IE angle
    kin(i,16)=RHF(2,1); %ML translation
    kin(i,17)=RHF(3,1)*cos(kin(i,13))-RHF(4,1)*sin(kin(i,13)); %AP translation (text has TF(1,4))
    kin(i,18)=(RHF(2,4)*RHF(2,1)+RHF(3,4)*RHF(3,1)+RHF(4,4)*RHF(4,1)); %SI translation
    
    p1 = gs(i,47:49); p2 = inter_rhip1(i,:);
    dist = sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2 + (p1(3) - p2(3))^2);
    deltaML = sin(kin(i,14))*dist;
    kin(i,16) = kin(i,16) + deltaML;
    
    
    %Right Ankle kinematics
    kin(i,19)=1*atan(RTA(3,4)/RTA(4,4)); %FE angle
    %        if kin(trial,19)<-.5 %To flip sign when trialon goes past 90 degrees
    %            kin(trial,19)=kin(trial,19)+pi;
    %        end
    kin(i,20)=acos(RTA(2,4))-pi/2; %VV angle
    kin(i,21)=atan(RTA(2,3)/RTA(2,2)); %IE angle
    kin(i,22)=RTA(2,1); %ML translation
    kin(i,23)=RTA(3,1)*cos(kin(i,19))-RTA(4,1)*sin(kin(i,19)); %AP translation (text has TF(1,4))
    kin(i,24)=(RTA(2,4)*RTA(2,1)+RTA(3,4)*RTA(3,1)+RTA(4,4)*RTA(4,1)); %SI translation
    
    p1 = gs(i,20:22); p2 = inter_rank1(i,:);
    dist = sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2 + (p1(3) - p2(3))^2);
    deltaML = sin(kin(i,20))*dist;
    kin(i,22) = kin(i,22) + deltaML;
        
end
%%
disp('test14')
%Write Kinematic Amplitude Card

f = kin; %kinematics
t = linspace(initial_time,final_time, size(f,1));
t = (t-initial_time); % set to start at 0s%=0.025/2
% t = initial_time:time_increment:final_time;
% t = (t-initial_time)*3; % set to start at 0s%=0.025/2

% f = kin;
% t = 0:time_increment:final_time;
% t = t*3;

fout1=fopen([path patient '\LtModel\amp_profiles_' name '.inp'],'w'); %

amp_list = {'AMP_RTF_FE_LOAD2'; 'AMP_RTF_VV_LOAD2'; 'AMP_RTF_IE_LOAD2'; 'AMP_RTF_ML_LOAD2'; 'AMP_RTF_AP_LOAD2'; 'AMP_RTF_SI_LOAD2'; ...
    'AMP_RPF_FE_LOAD2'; 'AMP_RPF_VV_LOAD2'; 'AMP_RPF_IE_LOAD2'; 'AMP_RPF_ML_LOAD2'; 'AMP_RPF_AP_LOAD2'; 'AMP_RPF_SI_LOAD2'; ...
    'AMP_RHip_FE_LOAD'; 'AMP_RHip_VV_LOAD'; 'AMP_RHip_IE_LOAD'; 'AMP_RHip_ML_LOAD'; 'AMP_RHip_AP_LOAD'; 'AMP_RHip_SI_LOAD'; ...
    'AMP_RAnk_FE_LOAD'; 'AMP_RAnk_VV_LOAD'; 'AMP_RAnk_IE_LOAD'; 'AMP_RAnk_ML_LOAD'; 'AMP_RAnk_AP_LOAD'; 'AMP_RAnk_SI_LOAD'};

for var = 1:size(amp_list,1)
    NAME = amp_list{var};
    fprintf(fout1,['*AMPLITUDE, VALUE=RELATIVE, SMOOTH=0.1, NAME=' NAME '\n']);
    for pc=1:size(f,1)
        if rem(pc,4)==1; fprintf(fout1,'%-1.4f,%-1.4f',[t(pc) f(pc,var)]);
        elseif rem(pc,4)==0; fprintf(fout1,',%-1.4f,%-1.4f\n',[t(pc) f(pc,var)]);
        else; fprintf(fout1,',%-1.4f,%-1.4f',[t(pc) f(pc,var)]);
        end
    end
    if rem(pc,4)>0; fprintf(fout1,'\n**\n'); %rem =remainder
    else
        fprintf(fout1,'**\n');
    end
end

fclose(fout1);
%%
%%
disp('test18')
%writes abaqus input file and output file replacing patient name, final time and output
%time increments

file=[path patient '/generic.inp'];
fin_r = fopen(file,'r');
TotTime=(final_time-initial_time); %from SLOW amp cards
interval_time=TotTime/100;
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
    else
        data{r}=temp;
        r=r+1;
    end
end

fin_w1 = fopen([path patient '\' name '.inp'],'w+');
for i=1:r-1
    fprintf(fin_w1,'%s',data{i});
    fprintf(fin_w1,'\n');
end

fclose(fin_w1);
