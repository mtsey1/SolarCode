function simedge(corr,edges,data,meta,k)
close all;
year=meta.Year;
res=1000;
days=365+~mod(year,4);
mesh=zeros(days, res);
sunarray=generate_sun_array(year,res,meta);
%might need to shift 0deg to south
% ----------------------------------
%add auto fit (only vert and horizontal to start
%==================================
sunarray(:,:,2)=mod(sunarray(:,:,2)+180,360);
leftheight=70;
rightheight=70;
leftedge=110;
rightedge=190;
centerheight=50;

for i=1:days
    for j=1:res
        if sunarray(i,j,3)>90 
            mesh(i,j)=1;
        elseif sunarray(i,j,2)<=leftedge && sunarray(i,j,3)>leftheight
            mesh(i,j)=1;
        elseif sunarray(i,j,2)>=rightedge && sunarray(i,j,3)>rightheight
            mesh(i,j)=1;
        elseif (sunarray(i,j,2)<rightedge && sunarray(i,j,2)>leftedge) && sunarray(i,j,3)>centerheight 
            mesh(i,j)=1;
        end 
    end
end
sunvec2=generate_sun_array(year,48,meta);
zenith2=sunvec2(:,11:41,3)';
azim2=mod(sunvec2(:,11:41,2)'+180,360);
zenith=sunarray(:,150:800,3)';
azim=sunarray(:,150:800,2)';
edge=diff(mesh');
[col,row]=find(edge);
col=col/1000*48-9;
X=[60, leftedge, leftedge, rightedge, rightedge,300];
Y=[leftheight,leftheight,centerheight,centerheight,rightheight,rightheight];
figure(100);
imagesc(squeeze(data(k,:,10:40))'); hold on;
plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
figure(101);
imagesc(squeeze(corr(k,:,:))'); hold on;
plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
figure(102);
imagesc(squeeze(edges(k,:,:,1)|edges(k,:,:,2))');hold on;
plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
axis([1 days 1 31])        
mornX=zenith2(squeeze(edges(k,:,:,1)|edges(k,:,:,2))');
mornY=azim2(squeeze(edges(k,:,:,1)|edges(k,:,:,2))');
figure(104);
plot(mornY,mornX,'*');
hold on
plot(azim(1:end,170),zenith(1:end,170));
hold on
plot(azim(1:end,350),zenith(1:end,350));
hold on
plot(X,Y);
axis([60 300 0 100])
set(gca, 'YDir','reverse');
figure
imagesc(mesh');
%itterative refining of paramaters
set(gca, 'YDir','reverse')
end
