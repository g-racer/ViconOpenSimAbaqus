function [time ] = getstoptimes(path, patient,trial, initial_time, final_time)
%function to find inital time when hit force plate to foot hitting ground
%and final time for full gait cycle based on knee angle. 
%have run iknverse kinematics for full cycle
force_file=[path patient '/ForceData/' trial '.mot'];
kin_file=[path patient '/IKresults/' trial '_ik.mot'];
%read in force data
fin_r = fopen(force_file,'r');
temp = fgetl(fin_r);
while strncmp(temp,'time',4) ~= 1
    temp = fgetl(fin_r);
end
p=1; n=1;
while ~feof(fin_r) %Do until the end of the file
    temp = fgetl(fin_r);
    if ~isempty(temp)
        ForceData(p,:) = str2num(temp);
        p=p+1;
    end
end

%read in ik data
fin_r = fopen(kin_file,'r');
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

time=struct;
 T=ForceData(:,1);
    Ktime=KinData(:,1); %different sampling freq crop force time and data
   [V, idxend] = min(abs(T - final_time));
    [V, idxstrt] = min(abs(T - initial_time));
    ind=idxstrt:idxend;
T=T(ind);

if  strncmp(trial,'lunge',5)==1
kneedata=KinData(:,13); % knee angle
[pks,locs_on] = min(kneedata); %find max knee flexion- angle reversed in opensim
time.maxflex=Ktime(locs_on);
end

% [pks,locs_on] = max(kneedata);
% time.minflex=ktime(locs_on);



if contains(trial,'gait')==1
    %FP 1 2 order =hits fp 2 then fp 1. 
   time.heelstrike=initial_time+(final_time-initial_time)*.05 ;
    forcedata=[ForceData(ind,3) ForceData(ind,10) ForceData(ind,17)]; %FP 2 and 3
   forcedata_s=smoothdata(forcedata, 'gaussian',30); %gaussian filter
          %find peaks from FP 2 for initial HS
            [maxtab, mintab]=peakdet(forcedata_s(:,2), 70, T) ;
    time.peak1=maxtab(1,1);
     time.peak2=maxtab(2,1);
     [slopex]=gradient(forcedata_s(:,2));
       locs_off=find(slopex < -0.4); %finds when force data decreases rapidly
                 time.toeoff=T(locs_off(end))-(final_time-initial_time)*.05;
  
else
 forcedata=[ForceData(ind,3)]; %FP 1
   forcedata_s=smoothdata(forcedata, 'gaussian',20);  
    time.heelstrike=initial_time+(final_time-initial_time)*.05 ;
    
if contains(trial,'stairup')==1
   

            [maxtab, mintab]=peakdet(forcedata_s, 10, T) ;
    time.peak1=maxtab(2,1);
     time.peak2=maxtab(end,1);
     
     [slopex]=gradient(forcedata_s);
       locs_off=find(slopex < -0.4); %finds when force data decreases rapidly
                 time.toeoff=T(locs_off(end));
      
elseif contains(trial,'stairdown')==1
   

            [maxtab, mintab]=peakdet(forcedata_s, 15, T) ;
    time.peak1=maxtab(2,1);
     time.peak2=maxtab(end,1);
     
     [slopex]=gradient(forcedata_s);
       locs_off=find(slopex < -0.4); %finds when force data decreases rapidly
                 time.toeoff=T(locs_off(end));
elseif contains(trial,'lunge')==1

          %find peaks from FP 1 
            [maxtab, mintab]=peakdet(forcedata_s, 10, T) ;
    time.peak1=maxtab(1,1);
     time.peak2=maxtab(end,1);
     
     [slopex]=gradient(forcedata_s);
       locs_off=find(slopex < -0.3); %finds when force data decreases rapidly
                 time.toeoff=T(locs_off(end));
                 
  elseif contains(trial,'stepup')==1
    %FP 1 2 order =hits fp 2 then fp 1. 
  
            [maxtab, mintab]=peakdet(forcedata_s, 5, T) ;
    time.peak1=maxtab(1,1);
     time.peak2=maxtab([end-1],1);
     
     [slopex]=gradient(forcedata_s);
       locs_off=find(slopex < -0.4); %finds when force data decreases rapidly
                 time.toeoff=T(locs_off(end));
end
%plot(T, forcedata_s)

%     time.startgait=initial_time+(final_time-initial_time)*.05 ;
%     time.secondgait=initial_time+(final_time-initial_time)*1/3;
%     time.thirdgait=initial_time+(final_time-initial_time)*2/3;
%         ind=find(ktime==initial_time); %initial time
%              locs_off=find(forcedata(ind:end,1) ==0);%finds when force data =0
%         time.offgait=ktime(locs_off(1)+ind); %inital time=location of first peak
    %finds time for first 5% 33% 66% of gait cycle. 
end



end





