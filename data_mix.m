datalist=[22,27,35];
noise=datamelb(datalist,:,:);
solarlist=1:100;
meta.postcode=meta.postcode(solarlist);
solardata=datasyd(solarlist,:,:);
simulateddata=zeros(length(datalist)*(1+length(solarlist)),366,48);
temp=[];
for i=1:length(datalist)
    for j=1:length(solarlist)
        simulateddata((i-1)+(j-1)*length(datalist)+1,1:365,:)=squeeze(noise(i,:,:)+solardata(j,1:365,:));
        simulateddata((i-1)+(j-1)*length(datalist)+1,366,:)=squeeze(noise(i,mod(round(60*rand)+336,365)+1,:)+solardata(j,366,:));
    end
end
for k=1:100
    temp=[temp,k,k,k];
end
state.postcode=state.postcode(temp);