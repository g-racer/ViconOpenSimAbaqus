function [ initial_time, final_time ] = gettrialtime(kin_F_file,force_file, patient,trial)
%function to find inital time when hit force plate to foot hitting ground
%and final time for full gait cycle based on knee angle. 
%have run iknverse kinematics for full cycle

%read in force data
fin_r = fopen(force_file,'r');
temp = fgetl(fin_r);
while strncmp(temp,'time',4) ~= 1
    temp = fgetl(fin_r);
end
p=1; 
while ~feof(fin_r) %Do until the end of the file
    temp = fgetl(fin_r);
    if ~isempty(temp)
        ForceData(p,:) = str2num(temp);
        p=p+1;
    end
end

%read in ik data
fin_r = fopen(kin_F_file,'r');
temp = fgetl(fin_r);
while strncmp(temp,'time',4) ~= 1
    temp = fgetl(fin_r);
end
p=1; n=1;
while ~feof(fin_r) %Do until the end of the file
    temp = fgetl(fin_r);
    if ~isempty(temp)
        KinData(p,:) = str2num(temp);
        p=p+1;
    end
end


if strncmp(patient,'VA',2)==1
    %%
if strncmp(trial,'RWalk',5)==1
    forcedata=ForceData(:,9); %FP 2 in y direction
 
    kneedata=KinData(:,11);
elseif strncmp(trial,'RCut',4)==1 || strncmp(trial,'RPivot',5)==1
    forcedata=ForceData(:,3); %FP 1 in y direction
  
    kneedata=KinData(:,20);
elseif strncmp(trial,'LCut',4)==1 || strncmp(trial,'LPivot',5)==1
    forcedata=ForceData(:,9); %FP 2 in y direction
   
    kneedata=KinData(:,20); %left knee angle
end
Ftime=ForceData(:,1);
[pks,locs_on] = findpeaks(forcedata); %find peak where strike force plate
initial_time=time(locs_on(1)); %inital time=location of first peak

time=KinData(:,1);
sf=length(Ktime)/length(time); %scale factor for time
ix=round(locs_on(1)*sf,0); %index


[pks,locs_on] = findpeaks(kneedata);
ix=find(pks>=60); %find when knee angle flexed
gait_time=Ktime(locs_on(ix(2)))-Ktime(locs_on(ix(1)));
final_time=initial_time+gait_time;

%if heel strike/full cycle not completed -use final time
if final_time>Ktime(end)
    final_time=Ktime(end); 
end




elseif strncmp(patient,'jw',2)==1
    
    time=ForceData(:,1);
    Ktime=KinData(:,1); %different sampling freq 
    idx = find(abs(time - Ktime(end))<=.005);
time=time(1:idx(end));
   
% 
%     if contains(trial,'down')==1 %manual
%         initial_time=3;
%         final_time=8.32;
%     else
    
        %initial time -fluro trials
        if strncmp(trial,'stair',5)==1 || strncmp(trial,'squat',5)==1 ...
                || strncmp(trial,'step',4)==1 || strncmp(trial,'chair',5)==1 || strncmp(trial,'lunge',5)==1
            
            forcedata=ForceData([1:length(time)],3); %FP 1 in y direction
            forcedata_s=smoothdata(forcedata, 'gaussian',30); %gaussian filter
            
            [slopex]=gradient(forcedata_s);
            %         jump = find(diff([forcedata_s(1); forcedata_s]) > 1, 1, 'first');
           if strncmp(trial,'lunge',5)==1 || contains(trial, 'down')==1
               locs_on=find(forcedata_s-forcedata_s(1)>20); %shallow increase in slope
           else
            locs_on=find(slopex > 0.5);%finds when force data increases rapidly
           end
           initial_time=time(locs_on(1)); %inital time=location of first peak
            
            %final time -fluro trials
            if strncmp(trial,'step',4)==1 || strncmp(trial,'lunge',5)==1%|| strncmp(trial,'chair',5)==1
               
                forcedata_s=smoothdata(forcedata, 'gaussian',100);
                
                [maxtab, mintab]=peakdet(forcedata_s, 100, time) ;
                final_time=mintab(1,1)+.1;
                
            elseif  strncmp(trial,'lunge',5)==1 
                forcedata_s=smoothdata(forcedata, 'gaussian',100);
                
                locs_off=find(min(forcedata_s)==forcedata_s);
                final_time=time(locs_off);
            else
                
                locs_off=find(slopex < -0.2); %finds when force data decreases rapidly
                final_time=time(locs_off(end));
            end
            
            % plot(time, forcedata_s)
            %times -gait trials
        elseif contains(trial,'gait')==1 || strncmp(trial,'bouncy',6)==1 ...
                || contains(trial,'crouch')==1 || strncmp(trial,'med',3)==1
            
            forcedata=[ForceData(1:length(time),10) ForceData(1:length(time),17)]; %FP 2 and 3
            %forcedata_s=smoothdata(forcedata,'gaussian',500); %gaussian filter
            [slopex1]=gradient(forcedata(:,1));
            locs_on1=find(slopex1 > 0.4);%finds when force data increases rapidly
            initial_time=time(locs_on1(1)); %inital time=location of first peak
            
            [slopex2]=gradient(forcedata(:,2));
            locs_on2=find(slopex2 > 0.5); %finds when force data increases rapidly
            final_time=time(locs_on2(1));
            
            if initial_time> final_time
                %if hit force plates in opposite order
                [initial_time, final_time] = swap(final_time, initial_time);
            end
            
        end
        
%     end
    fclose('all');
            
  % [pks,locs] = findpeaks(forcedata_s,'threshold', 10)'MinPeakHeight',10,'MinPeakDistance',10);
% smallPeakIndexes = pks < 300; %remove peaks when less than 300N
% pks(smallPeakIndexes) = [] ; %Reject Y value of peaks below this threshold
% locs(smallPeakIndexes) = [] ; 

end



