function  barplot(T, group1, group2, group1label, group2label, Title, ylab, Cimp)
Tstruct=T;
Bardata=struct;
tc=1;

Mean=[];
STD=[];

for idx1=1:length(group1)
    g1=group1{idx1};
    
    for idx2=1:length(group2)
        g2=group2{idx2};
        
        rows = strcmp(Tstruct.([group1label]),g1)==1 & strcmp(Tstruct.([group2label]),g2)==1;
        datatemp=table2array(Tstruct(rows,5));
        S = std(datatemp);
        M=mean(datatemp);
        Bardata.mean(tc)=M;
        Bardata.STD(tc)=S;
        Bardata.Load(tc)={g1};
        Bardata.Trial(tc)={g2};
        tc=tc+1;
        
        Mean(idx1,idx2)=M;
        STD(idx1,idx2)=S;
    end
    
    
end
impcolors={'b', 'y', 'r', 'g', 'm'};
impcolors= [0, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.8500, 0.3250, 0.0980;
    0.4660, 0.6740, 0.1880;
    0.4940, 0.1840, 0.5560	];
%
% 'b', 'y', 'r', 'g', 'm';'b', 'y', 'r', 'g', 'm';...
%     'b', 'y', 'r', 'g', 'm';'b', 'y', 'r', 'g', 'm';};

figure
%g2=implant
ax = axes;

h = bar(Mean,'BarWidth',1);

for k=1:5
    h(k).FaceColor = 'flat';
    
    c=repmat(impcolors(k,:),5,1);
    h(k).CData=c;
end

ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,1:length(group1));
% Naming each of the bar groups
xticklabels(ax,group1);
xlabel (group1label);
ylabel (ylab);
title(Title,'Interpreter', 'none')
% Creating a legend and placing it outside the bar plot
lg = legend(group2,'AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
ngroups = size(Mean, 1);
nbars = size(Mean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
X=zeros(ngroups, nbars);
%X=[];
for i = 1:nbars %for each implant
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Mean(:,i),STD(:,i), 'k', 'linestyle', 'none');
    X(:,i)=x'; %colums=each implant
end
X=reshape(X',25,1);
M=reshape(Mean',25,1);
st=reshape(STD',25,1);

if nargin ==8 %compare each group within implants
    Ci=Cimp([1:4 25:27 48:49 70 111:114 130:132 148:149 165 ...
        196:199 210:212 223:224 235 256:259 265:267 273:274 280 291:300],:);
    p=Ci(:,6); %Pvalues
      nos=sum(p<0.8);
      ind=find(p<0.8);
        xposi=zeros(1,nos);
          xposj=zeros(1,nos);
            yposi=zeros(1,nos);
            yposi=zeros(1,nos);
            P=zeros((2*nos));
             pct=1;
             
             Cip=Ci(ind,:);
              pval=Cip(:,6);
    for k=1:size(Cip,1)
    grpi=Cip(k,1);
       grpj=Cip(k,2);
            
            xposi(k)=X(grpi);
            xposj(k)=X(grpj);
            yposi(k)=M(grpi);
             yposj(k)=M(grpj) ; 
               erri(k)=st(grpi);
             errj(k)=st(grpj) ; 
             P(k,k+size(Cip,1))=pval(k);
  P(k+size(Cip,1),k)=pval(k);
    end
   P(P == 0) = NaN;
   X=[xposi xposj];
   Y=[yposi yposj];
   E=[erri errj];

             
            ypos=max([M(grpi) M(grpj)]) +min([M(grpi) M(grpj)]);
            pval(pct)=Ci(k,6);
            mysigstar(ax, [xposi xposj], ypos, pval);
            %
        end
        %
        %
    end