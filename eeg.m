data=readtable("S1_21_Male.csv");
data=table2array(data);
photoVector=data(:,17);
t=linspace(0,100,length(photoVector));
% figure (1); plot(t,photoVector);

th=240;
targetVector=double((photoVector>th));
 figure(2); plot(t,targetVector);

[pks,locs,w,p] = findpeaks(targetVector);

targetEEGs={};
channel16Vector=data(:,7);
for i=1:length(locs)
    targetEEGs{i}=channel16Vector(locs(i):locs(i)+w(i)).';
end

for i=1:length(targetEEGs)
    targetEEGs{i}=[targetEEGs{i},nan(1,max(w)-length(targetEEGs{i}))];
end

b=cell2mat(targetEEGs);
c=mean(b,1,'omitnan');

d=linspace(0,100,length(c))
figure(3); plot(d,c);
