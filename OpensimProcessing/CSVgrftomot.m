function [  ] = CSVgrftomot( patient, path, trial )
%Converts .csv grf data to Opensim .mot file with appropriate format.
%Based on Vicon to Opensim code.
name =[patient '_' trial];

forcedata_folder= [path patient '/ForceData/'];

csv_file = [forcedata_folder name '_grf.csv'] ;
mot_file = [forcedata_folder trial '.mot'] ;
%trc_file = [csv_folder name 'new.trc'] ;
signal = 0; no_fp = 0; no_mc = 1;


%
fin_r = fopen(csv_file,'r');
temp = fgetl(fin_r);
csv.head_grf{1} = temp;
p = 1; n = 1;
while ~feof(fin_r) %Do until the end of the file
    temp = fgetl(fin_r);
    k = strfind(temp,',,');
    % if ~isempty(k); temp = temp(1:k(1)-1); end
    if ~isempty(k) && signal == 0; temp = temp(1:k(1)-1);
    end
    if ~isempty(temp)
        
        csv.data_grf_str{p} = temp;
        csv.data_grf_num{1}(p,:) = str2num(temp);
        p=p+1;
    
        
    end
end
fclose(fin_r);


n_fp = (size(csv.data_grf_num{1},2)-1)/7;

cop_osim_col = cell(1,n_fp);
force_osim_col = cell(1,n_fp);
moment_osim_col = cell(1,n_fp);
for i = 1:3
    %COP
    temp = csv.data_grf_num{1}(:,i*7-2:(i*7));
    cop_osim_col{i}(:,3)=temp(:,1); %vicon x becomes opensim z
    cop_osim_col{i}(:,1)=temp(:,2); %vicon y becomes opensim x
    cop_osim_col{i}(:,2)=temp(:,3); %vicon z becomes opensim y
    %force
    temp = csv.data_grf_num{1}(:,(i*7-5):(i*7-3));
    force_osim_col{i}(:,3)=temp(:,1); %vicon x becomes opensim z
    force_osim_col{i}(:,1)=temp(:,2); %vicon y becomes opensim x
    force_osim_col{i}(:,2)=temp(:,3); %vicon z becomes opensim y
    
    moment_osim_col{i}= csv.data_grf_num{1}(:,(i*7+1));
    
    
    %         force_osim_col{i} = csv.data_grf_num{1}(:,(i*7-5):(i*7-3));
    %         moment_osim_col{i} = csv.data_grf_num{1}(:,(i*7+1));
end
time_grf_col = csv.data_grf_num{1}(:,1);

fin_w1 = fopen(mot_file,'w+');

%fprintf(fin_w1,'Ground reaction force (.MOT) by Navacchia Alessandro\n');
fprintf(fin_w1,[trial '.mot\n']);
fprintf(fin_w1,'version=1\n');

%     fprintf(fin_w1,'nRows=%d\n',size(csv.data_grf_num{1},1));
%     fprintf(fin_w1,'nColumns=%d\n',size(csv.data_grf_num{1},2)-1);
%        fprintf(fin_w1,'inDegrees=yes\n');

%     fprintf(fin_w1,'\n');
%   fprintf(fin_w1,['name' csv_file_name '.mot' '\n']);
% fprintf(fin_w1,'nRows=%d\n',size(csv.data_grf_num{1},1));
%     fprintf(fin_w1,'nColumns=%d\n',size(csv.data_grf_num{1},2)-1);
%   fprintf(fin_w1,'\n');
fprintf(fin_w1,'datacolumns %d\n',size(csv.data_grf_num{1},2));
fprintf(fin_w1,'datarows %d\n',size(csv.data_grf_num{1},1));
fprintf(fin_w1,'range %f %f\n',time_grf_col(1),time_grf_col(end));
fprintf(fin_w1,'endheader\n');
fprintf(fin_w1,'time\t');
for i = 1:n_fp
    fprintf(fin_w1,'%d_ground_force_vx\t',i);
    fprintf(fin_w1,'%d_ground_force_vy\t',i);
    fprintf(fin_w1,'%d_ground_force_vz\t',i);
    fprintf(fin_w1,'%d_ground_force_px\t',i);
    fprintf(fin_w1,'%d_ground_force_py\t',i);
    fprintf(fin_w1,'%d_ground_force_pz\t',i);
     fprintf(fin_w1,'%d_ground_torque_y\t',i);
end
% for i = 1:n_fp
%     fprintf(fin_w1,'%d_ground_torque_y\t',i);
% end
fprintf(fin_w1,'\n');

for j = 1:size(time_grf_col,1)
    fprintf(fin_w1,'%f\t',time_grf_col(j));
    for i = 1:n_fp
        fprintf(fin_w1,'%f\t',force_osim_col{i}(j,1));
        fprintf(fin_w1,'%f\t',force_osim_col{i}(j,2));
        fprintf(fin_w1,'%f\t',force_osim_col{i}(j,3));
        fprintf(fin_w1,'%f\t',cop_osim_col{i}(j,1));
        fprintf(fin_w1,'%f\t',cop_osim_col{i}(j,2));
        fprintf(fin_w1,'%f\t',cop_osim_col{i}(j,3));
        fprintf(fin_w1,'%f\t',moment_osim_col{i}(j,1));
    end
    %         for i = 1:n_fp
    %             fprintf(fin_w1,'%f\t',moment_osim_col{i}(j,1));
    %             if isnan(moment_osim_col{i}(j,2))
    %                 fprintf(fin_w1,'%d\t',0);
    %             else fprintf(fin_w1,'%f\t',moment_osim_col{i}(j,2));
    %             end
    %             fprintf(fin_w1,'%f\t',moment_osim_col{i}(j,3));
    %         end
    fprintf(fin_w1,'\n');
end
fclose('all');
disp(['GRF file for ' trial ' completed'])
end

