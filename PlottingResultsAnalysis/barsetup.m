function [X,Y,E,P,grps,Pr,Prgrps] =  barsetup(T, gr1_Trial, gp2_Implant, group1label, group2label, Cimp, grp)
Tstruct=T;
Bardata=struct;
tc=1;

Mean=[];
STD=[];

for idx1=1:length(gr1_Trial)
    g1=gr1_Trial{idx1};
    
    for idx2=1:length(gp2_Implant)
        g2=gp2_Implant{idx2};
        
        rows = strcmp(Tstruct.([group1label]),g1)==1 & strcmp(Tstruct.([group2label]),g2)==1;
        datatemp=table2array(Tstruct(rows,5));
        S = std(datatemp);
        M=mean(datatemp);
        %coloumn =imp
        Mean(idx1,idx2)=M;
        STD(idx1,idx2)=S;
    end
    
    
end
Y=Mean;
E=STD;

ngroups = size(Mean, 1);
nbars = size(Mean, 2);
% Calculating the width for each bar group

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
X=zeros(ngroups, nbars);
groupwidth = min(0.8, nbars/(nbars + 1.5));

%X=[];
for i = 1:nbars %for each implant
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    X(:,i)=x; %colums=each implant
end

if nargin ==7 %compare each group within implants
    P=zeros(size(grp,1),size(grp,1));
    Pr=zeros(10,6);
    grps=cell(size(grp,1),size(grp,1));
    Prgrps=cell(10,6);
    ol=1;
    for k=1:length(Cimp(:,1))
       
        if Cimp(k,6)<=.05
             row=Cimp(k,1);
        col=Cimp(k,2);%get group numbers
            if strcmp( grp{col}(end-7:end),grp{row}(end-7:end))==1 %compare if from same trial
                trial=grp{col}(end-7:end);
                Cimpfill(ol,:)=Cimp(k,[1:2 6]);
          
                ol=ol+1;
                P(row,col)=Cimp(k,6);
                P(col,row)=Cimp(k,6);
                grps{col,row}=([grp{col} grp{row}]);
                ind5=5:6:30;
                ind6=6:6:30;
                
             p5=P(ind5,ind5) ;  
              p6=P(ind6,ind6) ;  
                switch trial %row no.
                    case '=lungef1'; r=1;
                    case 'gait_og1'; r=2;
                    case 'irdownf2'; r=3;
                    case  'tairupf1'; r=4;
                    case  'stepupf2'; r=5;
                end
                if strcmp((grp{col}(12)),'5')==1 && strcmp((grp{row}(12)),'6')==1 ...
                        || strcmp((grp{col}(12)),'6')==1 && strcmp((grp{row}(12)),'5')==1
                    %if implant 5 & 6 different
                    imp=min([str2num(grp{col}(12)) str2num(grp{row}(12))]); %implant no. 
                    Pr(r+5,imp)=Cimp(k,6);
                    Prgrps{r+5,imp}=([grp{col} grp{row}]);
                    
                elseif strcmp((grp{col}(12)),'5')==1 || strcmp((grp{row}(12)),'5')==1
                    imp=min([str2num(grp{col}(12)) str2num(grp{row}(12))]);
                    Pr(r,imp)=Cimp(k,6);
                    Prgrps{r,imp}=([grp{col} grp{row}]);
                    
                elseif strcmp((grp{col}(12)),'6')==1 || strcmp((grp{row}(12)),'6')==1 %if different from patient
                    imp=min([str2num(grp{col}(12)) str2num(grp{row}(12))]);
                    Pr(r+5,imp)=Cimp(k,6);
                    Prgrps{r+5,imp}=([grp{col} grp{row}]);
                end
                
            end
        end
    end
    P(P == 0) = NaN;
    Pr(Pr == 0) = NaN;
end



%     Ci=Cimp([1:4 25:27 48:49 70 111:114 130:132 148:149 165 ...
%         196:199 210:212 223:224 235 256:259 265:267 273:274 280 291:300],:);
%     p=Ci(:,6); %Pvalues
%       nos=sum(p<0.8);
%       ind=find(p<0.8);
%
%             P=zeros(25);
%              pct=1;
%
%              Cip=Ci(ind,:);
%               pval=Cip(:,6);
%     for k=1:size(Cip,1)
%     grpi=Cip(k,1);
%        grpj=Cip(k,2);
%
%
%              P(grpi,grpj)=pval(k);
%   P(grpj,grpi)=pval(k);
%     end
%    P(P == 0) = NaN;
% end
% end

%
%

end
