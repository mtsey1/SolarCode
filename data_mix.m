datalist=[4,7,13,14,17,19,20,21,22,23,25,26,27,29,31,35];
noise=data(datalist,:,:);
solarlist=[6,8,9,10,13,14,15,19,20,21,24,26,27,30,35];
solardata=datasyd(solarlist,:,:);
simulateddata=zeros(length(datalist)*(1+length(solarlist)),366,48);
for i=1:length(datalist)
    for j=1:length(solarlist)
        simulateddata((i-1)+(j-1)*length(datalist)+1,1:365,:)=squeeze(noise(i,:,:)+solardata(j,1:365,:));
        %simulateddata(i+j*length(datalist),366,:)=data(i,mod(round(60*rand)+336,366),:)+datasyd(j,366,:);
    end
end
