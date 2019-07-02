function res=manipulate(res,correlation,data,meta,edges);

year=meta.Year;
lengthtot=100;
plotfind=0;
listin=0;
%logistic regression

%method to identify house from plot
%(potentially scatter3 with particular view)
if listin 
    houses=listHouse;
else
    houses=1:lengthtot;
end
sunvec=generate_sun_array(year,48,meta);
mornshaded=res(res(houses,9,1)==1);
aftershaded=res(res(houses,9,2)==1);
mornnonshaded=res(res(houses,9,1)==0);
afternonshaded=res(res(houses,9,2)==0);
figure(1); plot(res(mornshaded,10,1),res(mornshaded,11,1),'ro'); hold on;
plot(res(aftershaded,10,2),res(aftershaded,11,2),'ro'); hold on;
plot(res(mornnonshaded,10,1),res(mornnonshaded,11,1),'b*'); hold on;
plot(res(afternonshaded,10,2),res(afternonshaded,11,2),'b*');
xlabel('edge angle comparison')
ylabel('deviation within region')

zenith=sunvec(:,11:41,3)';
azim=mod(sunvec(:,11:41,2)'+180,360);

figure(2); plot(res(mornshaded,13,1),res(mornshaded,11,1),'ro'); hold on;
plot(res(aftershaded,13,2),res(aftershaded,11,2),'ro'); hold on;
plot(res(mornnonshaded,13,1),res(mornnonshaded,11,1),'b*'); hold on;
plot(res(afternonshaded,13,2),res(afternonshaded,11,2),'b*');
xlabel('edge strength')
ylabel('deviation within region')
a=[res(mornshaded,1,1),res(aftershaded,1,1),res(mornnonshaded,1,1),res(aftershaded,1,1)];
if plotfind
    [xx,yy,ind]=graphpoints;
    houses=a(ind)
end
    
for i=houses
    figure(100);
    imagesc(squeeze(data(i,:,:))');
    figure(101);
    imagesc(squeeze(correlation(i,:,:))')
    figure(102);
    imagesc(squeeze(edges(i,:,:,1)|edges(i,:,:,2)));
    axis([1 days 1 31])        
    mornX=zenith(squeeze(edges(i,:,:,1)|edges(i,:,:,2)));
    mornY=azim(squeeze(edges(i,:,:,1)|edges(i,:,:,2)));
    figure(104);
    plot(mornY,mornX,'*');
    hold on
    plot(azim(1:end,170),zenith(1:end,170));
    hold on
    plot(azim(1:end,350),zenith(1:end,350));
    axis([60 300 0 100])
    set(gca, 'YDir','reverse');
    if vischeck
        res(i,9,1)=input('Morning shading True False\n');
        res(i,9,2)=input('Afternoon shading True False\n');

    end
    update=input('Update plots and logistic regression\n');
    if update
        close([1,2])
        fitdat=[res(1:lengthtot,[10,11,13],1),res(1:lengthtot,[10,11,13],2)];
        [B,dev,stats]=mnrfit(fitdat,[res(1:lengthtot,9,1),res(1:lengthtot,9,2)]);
        mornshaded=res(res(houses,9,1)==1);
        aftershaded=res(res(houses,9,2)==1);
        mornnonshaded=res(res(houses,9,1)==0);
        afternonshaded=res(res(houses,9,2)==0);
        figure(1); plot(res(mornshaded,10,1),res(mornshaded,11,1),'ro'); hold on;
        plot(res(aftershaded,10,2),res(aftershaded,11,2),'ro'); hold on;
        plot(res(mornnonshaded,10,1),res(mornnonshaded,11,1),'b*'); hold on;
        plot(res(afternonshaded,10,2),res(afternonshaded,11,2),'b*');
        xlabel('edge angle comparison')
        ylabel('deviation within region')

        figure(2); plot(res(mornshaded,13,1),res(mornshaded,11,1),'ro'); hold on;
        plot(res(aftershaded,13,2),res(aftershaded,11,2),'ro'); hold on;
        plot(res(mornnonshaded,13,1),res(mornnonshaded,11,1),'b*'); hold on;
        plot(res(afternonshaded,13,2),res(afternonshaded,11,2),'b*');
        xlabel('edge strength')
        ylabel('deviation within region')
    end
    fprintf(repmat('\b',1,72));
end
end
    
   