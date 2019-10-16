%display function
house=38;
close all;
s=state;
%data plot 
year=meta.Year;
daystot=365+~mod(year,4);
[sun,panel]=sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(house),s.solar_ze(house),0,s.dark_end,meta.Days);  
sunvec=generate_sun_array(year,48,meta);
zenith2=squeeze(sunvec(:,11:41,3)');
azim2=squeeze(mod(sunvec(:,11:41,2)'+180,360));
row=simed(house,:,1);
col=simed(house,:,2);
mornedg=squeeze(edges(house,:,:,1));
afteredg=squeeze(edges(house,:,:,2));
days=1:daystot;
figure(100);
    plot(days, panel(1,:)-10, days, panel(2,:)-10, 'LineWidth', 2)
    title('Sunrise and Sunset')
    xlabel('Day of Year')
    ylabel('Time')
imagesc(squeeze(data(house,:,10:40))'); hold on;
%seperate top and bottom curve        
plot(row(col<=12),col(col<=12),row(col>17),col(col>17));
figure(101);
imagesc(squeeze(correlate(house,:,:))'); hold on;
plot(row(col<=12),col(col<=12),row(col>17),col(col>17));
figure(102);
imagesc((mornedg|afteredg)');hold on;
plot(row(col<=12),col(col<=12),row(col>17),col(col>17));
%axis([1 days 1 31])
mornX=zenith2((mornedg|afteredg)');
mornY=azim2((mornedg|afteredg)');
figure(104);
plot(mornY,mornX,'*');
hold on
plot(azim2(1:end,170),zenith2(1:end,170));
hold on
plot(azim2(1:end,350),zenith2(1:end,350));
hold on
xmorn=180:300;
xafter=60:180;
ymorn=fitvals(house,1,1)*xmorn+fitvals(house,1,2);
yafter=fitvals(house,2,1)*xmorn+fitvals(house,2,2);
plot(xmorn,ymorn,xafter,yafter)
axis([60 300 0 100])
set(gca, 'YDir','reverse');

