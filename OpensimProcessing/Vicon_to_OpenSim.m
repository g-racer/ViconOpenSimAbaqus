% Function written by Alessandro Navacchia, University of Denver

% variables to try the code out
%csv_file = 'P004 pre raw.csv'
%mot_file = 'CTR_pilot.mot'
%trc_file = 'CTR_pilot.trc'
%direction = 1
%vicon_version = 'CTR'
% 
% The inputs needed by this function are:
% - the csv_file location (obtained with Vicon Nexus extracting both the GRF and the markers locations)
% - the output files locations (trc_file for the markers and mot_file for the GRF data)
% - the Vicon Nexus version (that is the number x if the version is 1.8.x) this code currently works with 1.8.2. Vicon version 2.3 for COBR and CTR for CU data. 
% - the direction of motion (1 if you want the OpenSim positive anterior direction to be oriented 
%       as the Vicon positive anterior direction, 0 if the opposite)

function Vicon_to_OpenSim(csv_file,trc_file,mot_file,vicon_version,direction)
%% Read csv_file and save headers and data in structure csv
tic 

if direction == 1
    sign_dir = 1;
else sign_dir = -1;
end


[~,csv_file_name,~] = fileparts(csv_file);
%
fin_r = fopen(csv_file,'r');
temp = fgetl(fin_r);
if strncmp(temp,'Trajec',6) == 1; signal = 1; no_fp = 1; no_mc = 0;
else signal = 0; no_fp = 0; 
end

if no_fp == 0
    
    while strncmp(temp,'Devices',7) ~= 1
        temp = fgetl(fin_r);
    end
    
    csv.head_grf{1} = temp;
    for row = 2:5
        csv.head_grf{row} = fgetl(fin_r);
    end
    p = 1; n = 1;
    while ~feof(fin_r) %Do until the end of the file
        temp = fgetl(fin_r);
        k = strfind(temp,',,');
        % if ~isempty(k); temp = temp(1:k(1)-1); end
        if ~isempty(k) && signal == 0; temp = temp(1:k(1)-1); end
        if ~isempty(k) && signal == 1 
            for s = 1:length(k)
                temp = [temp(1:k(s)) 'NaN' temp(k(s)+1:end)];
                if s ~= length(k); k(s+1:end) = k(s+1:end)+3; end
            end
        end
        if ~isempty(temp)
            if signal == 0
                csv.data_grf_str{p} = temp;
                csv.data_grf_num{1}(p,:) = str2num(temp);
                p=p+1;
            else
                if no_mc == 0
                    csv.data_mocap_str{n} = temp; % read in the next line
                    csv.data_mocap_num{1}(n,:) = str2num(temp);
                    n=n+1;
                end
            end
        else
            if signal == 0
                temp = fgetl(fin_r);
                
                while strncmp(temp,'Trajec',6) ~= 1
                    temp = fgetl(fin_r);
                end
                
                if strncmp(temp,'Trajec',6) == 1
                    csv.head_mocap{1} = temp;
                    for row = 2:5
                        csv.head_mocap{row} = fgetl(fin_r);
                    end
                    no_mc = 0;
                else no_mc = 1;
                end
            end
            signal = signal+1;
        end
    end
    fclose(fin_r);
else
    csv.head_mocap{1} = temp;
    for row = 2:5
        csv.head_mocap{row} = fgetl(fin_r);
    end
    fclose(fin_r);
    csv.data_mocap_num{1} = dlmread(csv_file,',',5,0);
end
if isfield(csv,'data_mocap_num') == 0
    no_mc = 1;
end
%
%% Calculate transformed GRF data for OpenSim coord system
if no_fp == 0
    if strcmp(vicon_version,'CTR')
% get rid of EMG columns
title_line = textscan(csv.head_grf{4},repmat('%s',[1,n]),'Delimiter',',');
j = 0;
for i = 1:length(title_line)
    if strncmp(title_line{i},'EMG',3) == 1
        break
    else j = j + 1;
    end
end
    
csv.data_grf_num{1} = csv.data_grf_num{1}(:,1:j);
    end
n_fp = (size(csv.data_grf_num{1},2)-2)/9;

cop_vicon_col = cell(1,n_fp);
force_vicon_col = cell(1,n_fp);
moment_vicon_col = cell(1,n_fp);

% frame_grf_col = csv.data_grf_num{1}(:,1);
frame_grf_col = csv.data_grf_num{1}(:,1) - min(csv.data_grf_num{1}(:,1)) + 1;
if vicon_version == 2
    for i = 1:n_fp
        cop_vicon_col{i} = csv.data_grf_num{1}(:,i*3:(i*3+2));
        % negative because when export from vicon signs are flipped!
        force_vicon_col{i} = -csv.data_grf_num{1}(:,(i*3+n_fp*3):(i*3+2+n_fp*3));
        moment_vicon_col{i} = -csv.data_grf_num{1}(:,(i*3+n_fp*6):(i*3+2+n_fp*6));
    end
elseif vicon_version == 5
    for i = 1:n_fp
        force_vicon_col{i} = csv.data_grf_num{1}(:,(9*i-6):(9*i-6+2));
        moment_vicon_col{i} = csv.data_grf_num{1}(:,(9*i-3):(9*i-3+2));
        cop_vicon_col{i} = csv.data_grf_num{1}(:,(9*i):(9*i+2));
    end
elseif strcmp(vicon_version,'2.3')
    for i = 1:n_fp
        force_vicon_col{i} = -csv.data_grf_num{1}(:,(9*i-6):(9*i-6+2));
        moment_vicon_col{i} = -csv.data_grf_num{1}(:,(9*i-3):(9*i-3+2));
        cop_vicon_col{i} = csv.data_grf_num{1}(:,(9*i):(9*i+2));
    end
elseif strcmp(vicon_version,'CTR')
    for i = 1:n_fp
        force_vicon_col{i} = -csv.data_grf_num{1}(:,(9*i-6):(9*i-6+2));
        moment_vicon_col{i} = -csv.data_grf_num{1}(:,(9*i-3):(9*i-3+2));
        cop_vicon_col{i} = csv.data_grf_num{1}(:,(9*i):(9*i+2));
    end
end

cop_osim_col = cell(1,n_fp);
force_osim_col = cell(1,n_fp);
moment_osim_col = cell(1,n_fp);

if strcmp(vicon_version,'CTR')
    for i = 1:n_fp
        force_osim_col{i} = [sign_dir*force_vicon_col{i}(:,1) force_vicon_col{i}(:,3) -sign_dir*force_vicon_col{i}(:,2)];
        moment_osim_col{i} = zeros(size(moment_vicon_col{i},1),3);
        h = moment_vicon_col{i}(1,1);
        xp = (h*force_vicon_col{i}(:,2) - moment_vicon_col{i}(:,1))./force_vicon_col{i}(:,3);
        yp = (-h*force_vicon_col{i}(:,1) - moment_vicon_col{i}(:,2))./force_vicon_col{i}(:,3);
        moment_osim_col{i}(:,2) = ((moment_vicon_col{i}(:,3) - xp.*force_vicon_col{i}(:,1) - yp.*force_vicon_col{i}(:,2)))/1000;
        cop_osim_col{i} = [sign_dir*cop_vicon_col{i}(:,1) cop_vicon_col{i}(:,3) -sign_dir*cop_vicon_col{i}(:,2)]/1000;
    end
else
    for i = 1:n_fp
        force_osim_col{i} = [sign_dir*force_vicon_col{i}(:,2) force_vicon_col{i}(:,3) sign_dir*force_vicon_col{i}(:,1)];
        moment_osim_col{i} = zeros(size(moment_vicon_col{i},1),3);
        h = moment_vicon_col{i}(1,2);
        xp = (-h*force_vicon_col{i}(:,1) - moment_vicon_col{i}(:,2))./force_vicon_col{i}(:,3);
        yp = (-h*force_vicon_col{i}(:,2) + moment_vicon_col{i}(:,1))./force_vicon_col{i}(:,3);
        moment_osim_col{i}(:,2) = ((moment_vicon_col{i}(:,3) - xp.*force_vicon_col{i}(:,2) + yp.*force_vicon_col{i}(:,1)))/1000;
        cop_osim_col{i} = [sign_dir*cop_vicon_col{i}(:,2) cop_vicon_col{i}(:,3) sign_dir*cop_vicon_col{i}(:,1)]/1000;
    end
end

end
%
%% Write transformed GRF data in mot_file
if no_fp == 0
    
    % find mocap freq
    line = textscan(csv.head_mocap{2},repmat('%d',[1,n]),'Delimiter',',');
    mocap_freq = double(line{1});
    % find grf freq
    line = textscan(csv.head_grf{2},repmat('%d',[1,n]),'Delimiter',',');
    grf_freq = double(line{1});
    % ratio between freq
    % freq_ratio = grf_freq/mocap_freq;
    
    fin_w1 = fopen(mot_file,'w+');
    
    fprintf(fin_w1,'Ground reaction force (.MOT) by Navacchia Alessandro\n');
    fprintf(fin_w1,'nRows=%d\n',size(csv.data_grf_num{1},1));
    fprintf(fin_w1,'nColumns=%d\n',size(csv.data_grf_num{1},2)-1);
    fprintf(fin_w1,'\n');
    fprintf(fin_w1,['name' csv_file_name '_grf.mot' '\n']);
    fprintf(fin_w1,'datacolumns %d\n',size(csv.data_grf_num{1},2)-1);
    fprintf(fin_w1,'datarows %d\n',size(csv.data_grf_num{1},1));
    fprintf(fin_w1,'range %f %f\n',(frame_grf_col(1)-1)/mocap_freq,frame_grf_col(end)/mocap_freq-0.0005);
    fprintf(fin_w1,'endheader\n');
    fprintf(fin_w1,'time\t');
    for i = 1:n_fp
        fprintf(fin_w1,'%d_ground_force_vx\t',i);
        fprintf(fin_w1,'%d_ground_force_vy\t',i);
        fprintf(fin_w1,'%d_ground_force_vz\t',i);
        fprintf(fin_w1,'%d_ground_force_px\t',i);
        fprintf(fin_w1,'%d_ground_force_py\t',i);
        fprintf(fin_w1,'%d_ground_force_pz\t',i);
    end
    for i = 1:n_fp
        fprintf(fin_w1,'%d_ground_torque_x\t',i);
        fprintf(fin_w1,'%d_ground_torque_y\t',i);
        fprintf(fin_w1,'%d_ground_torque_z\t',i);
    end
    fprintf(fin_w1,'\n');
    
    frame_length = 1.0/grf_freq;
    for j = 1:size(frame_grf_col,1)
        fprintf(fin_w1,'%f\t',double((frame_grf_col(1)-1)/mocap_freq)+(j-1)*frame_length);
        for i = 1:n_fp
            fprintf(fin_w1,'%f\t',force_osim_col{i}(j,1));
            fprintf(fin_w1,'%f\t',force_osim_col{i}(j,2));
            fprintf(fin_w1,'%f\t',force_osim_col{i}(j,3));
            fprintf(fin_w1,'%f\t',cop_osim_col{i}(j,1));
            fprintf(fin_w1,'%f\t',cop_osim_col{i}(j,2));
            fprintf(fin_w1,'%f\t',cop_osim_col{i}(j,3));
        end
        for i = 1:n_fp
            fprintf(fin_w1,'%f\t',moment_osim_col{i}(j,1));
            if isnan(moment_osim_col{i}(j,2))
                fprintf(fin_w1,'%d\t',0);
            else fprintf(fin_w1,'%f\t',moment_osim_col{i}(j,2));
            end
            fprintf(fin_w1,'%f\t',moment_osim_col{i}(j,3));
        end
        fprintf(fin_w1,'\n');
    end
    fclose(fin_w1);
end
%
%% Calculate transformed Motion Capture data for OpenSim coord system
if no_mc == 0
n_data = sum(~isnan(csv.data_mocap_num{1}(1,:)))+sum(isnan(csv.data_mocap_num{1}(1,:))); 
csv.data_mocap_num{1} = csv.data_mocap_num{1}(:,[1 3:n_data]); %ignore second column
n_markers = (size(csv.data_mocap_num{1},2)-1)/3; 


% frame_mocap_col = csv.data_mocap_num{1}(:,1);
% subframe_mocap_col = csv.data_mocap_num{1}(:,2);

loc_vicon_col = cell(1,n_markers);
for i = 1:n_markers
    loc_vicon_col{i} = csv.data_mocap_num{1}(:,i*3-1:(i*3+1));
end

loc_osim_col = cell(1,n_markers);

if strcmp(vicon_version,'CTR')
    for i = 1:n_markers
        loc_osim_col{i} = [sign_dir*loc_vicon_col{i}(:,1) loc_vicon_col{i}(:,3) -sign_dir*loc_vicon_col{i}(:,2)];
    end
else
    for i = 1:n_markers
        loc_osim_col{i} = [sign_dir*loc_vicon_col{i}(:,2) loc_vicon_col{i}(:,3) sign_dir*loc_vicon_col{i}(:,1)];
    end
end
end
%
%% Write transformed Motion Capture data in trc_file
if no_mc == 0
temp = csv.head_mocap{2}(1:3);
k = strfind(temp,',,');
if ~isempty(k); temp = temp(1:k(1)-1); end
rate = str2num(temp);
% startFrame = csv.data_mocap_num{1}(1,1);
startFrame = 1;
numRows = size(csv.data_mocap_num{1},1);
markers = textscan(csv.head_mocap{3},'%s','delimiter',',');

fin_w2 = fopen(trc_file,'w+');

fprintf(fin_w2,['PathFileType\t4\t(X/Y/Z)\t' csv_file_name '.trc\tMaker trajectories (.TRC) Navacchia Alessandro\n']);
fprintf(fin_w2,'DateRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStarFrame\tOrigNumFrames\n');
fprintf(fin_w2,'%d\t%d\t%d\t%d\tmm\t%d\t%d\t%d\n',rate,rate,numRows,n_markers,rate,startFrame,numRows);
fprintf(fin_w2,'Frame#\tTime\t');
for i = 1:n_markers
    temp = markers{1}{i*3,1};
    index = strfind(temp,':');
    fprintf(fin_w2,temp(index+1:end));
    fprintf(fin_w2,'\t\t\t');
end
fprintf(fin_w2,'\n');
fprintf(fin_w2,'\t\t');
for i = 1:n_markers
    fprintf(fin_w2,'X%d\tY%d\tZ%d\t',i,i,i);
end
fprintf(fin_w2,'\n');
fprintf(fin_w2,'\n');

for j = 1:numRows
    fprintf(fin_w2,'%d\t',startFrame+j-1);
    fprintf(fin_w2,'%f\t',double(startFrame+j-2)/double(mocap_freq));
    for i = 1:n_markers
        if isnan(loc_osim_col{i}(j,1)) || loc_osim_col{i}(j,1) == 0; fprintf(fin_w2,'\t');
        else fprintf(fin_w2,'%f\t',loc_osim_col{i}(j,1));
        end
        if isnan(loc_osim_col{i}(j,2)) || loc_osim_col{i}(j,2) == 0; fprintf(fin_w2,'\t');
        else fprintf(fin_w2,'%f\t',loc_osim_col{i}(j,2));
        end
        if isnan(loc_osim_col{i}(j,3)) || loc_osim_col{i}(j,3) == 0; fprintf(fin_w2,'\t');
        else fprintf(fin_w2,'%f\t',loc_osim_col{i}(j,3));
        end
    end
    fprintf(fin_w2,'\n');
end

fclose(fin_w2);
end
toc
end
%
%