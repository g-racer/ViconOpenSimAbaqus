function [  ] = trcmarkerformat( patient, path, trial ) 
%formats the csv file to remove decimal places after the frame number and
%add in a newline after the headers so it is compatible with Opensim.
markerdata_folder = [path patient '/MarkerData/'];
trc_file=[path patient '/Synchronized Motion Data/Fluoroscopy Trials/Video Motion Data/' patient '_' trial '_new.trc'];
  trc_newfile=[markerdata_folder trial '.trc']
fid = fopen(trc_file,'r');
headers=cell(3,1);
for i=1:5
    headers{i} = fgetl(fid);
    %if contains(headers{i},'txt')
  headers{i}= replace(headers{i},'txt','trc  ');
  headers{i}= replace(headers{i},[patient '_' trial],trial);
    %end
end

trc.data_mocap_num{1} = dlmread(trc_file,'\t',5,0);
%trc.data_mocap_num{1}(1,:)=[];
trc.data_mocap_num{1}(:,end)=[]; %coloumn of 0 at end
%modifies coodinate system -comment if dont need. 
temp=dlmread(trc_file,'\t',5,0);
nmarkers=(size(trc.data_mocap_num{1},2)-3)/3;
for i=1:nmarkers
    trc.data_mocap_num{1}(:,3*i+2)=temp(:,3*i); %vicon x becomes opensim z
     trc.data_mocap_num{1}(:,3*i)=temp(:,3*i+1); %vicon y becomes opensim x
 trc.data_mocap_num{1}(:,3*i+1)=temp(:,3*i+2); %vicon z becomes opensim y
end
 nmrk=(size(trc.data_mocap_num{1},2)-3)/3;
%trc.data_mocap_num{1}(:,2:end)=linspace(1,pi/2,size(trc.data_mocap_num{1},1));
fid = fopen(trc_newfile,'w');
for i=1:5
    fprintf(fid,headers{i});
    fprintf(fid,'\n');
end
fprintf(fid,'\n'); %add in new line

for i = 1:size(trc.data_mocap_num{1},1)
    fprintf(fid, '%.0f\t',round(trc.data_mocap_num{1}(i,1))); %remove decimal places for frame no.
    fprintf(fid,'%f\t',trc.data_mocap_num{1}(i,[2:end]));
    fprintf(fid,'\n');
end
fclose('all');
disp([trial ' marker formatting completed'])
end

