function [correlation,res,edges,meta]=cloudcorr(s,data,meta) 
%------------------------------------------------------------
%cloudcorr takes in the gross power data for a household and the meta 
%and state data for the household it then conducts a rolling correlation
%between the data and the cloud cover (contained within s.solar_by_pc).
%INPUTS ----------------------------------------------------------
% s: state structure data generated from phase 1 and phase 2 data which includes the
% capacity, postcode, cloud cover and other data from the dataset
% data: this is the 3 dimentional array of household data for each house
% (dimension 1) for each day (dimension2) and for each half hour period
% (dimension 3)
% meta: contains all the meta data for the dataset including positional
% measurements, year.
%VARIABLES 
%----------------------------------------------------------------------
n=60; %length of correlation
width=3; % width of correlation
days=meta.Days; %number of days in the year
total=8000; %Max number of household analysed  
showplot=0; % logical value as to whether to show plots and figures generated 
shadowonly=0; %logical to use subset of data specified by vector shadowinglist
corrmanip=1; %logical to manipulate correlation plot around sunrise/set
noisered=1; %logical to reduce the amplitude of usage noise
vischeck=0; %logical value allowing to input the 
%B: old logistic regrestion values
B=[10.8481630139160;-16.2828877554926;-1.06561285029183;5.08033613741109];
correlationtrig=-0.35; %threshold 
%shadowinglist: a list of specific houses that are needed to be checked 
shadowinglist=[20,146,179,208,233,238,239,240,243,259,266,277,284,292];
%OUTPUTS----------------------------------------------------------------
%Correlation: array of correlation post procesing between the power data
%and the cloud cover data within the the 10th to 40th half hour segments of
%the day. it has dimensions of house number (dimension 1) day of the year 
%(dimension2) and half hour segments begining from the 10th half hour segment
%Result: is a 3 dimensional vector showing the quality of the shading edge 
%the first dimension is the house number. the second dimension is the list
%of different paramaters which are as follows
%   1: house number
%   2-8: are being used via the threshold calculations and isnt current prefered method   
%   2: ratio of days with shading to total days
%   3: smoothness of the edge 
%   4: edge angle comparison 
%   5: average of previous measurements
%   6: threshold of data
%   7: threshold of correlation 
%   8: **empty**
%   9: either manual shading edge true false or the calculated probability
%   from logistic regression
%   10-13: these values come from the quality metrics calculated from the canny edge detection
%   10: edge angle comparison
%   11: deveation of lowest 30% of data within shaded region
%   12,13: tresholds from the canny edge detection
%the 3rd dimension is morning (row 1) and afternoo row 2
%-------------------------------------------------------------------
%Edited 21/02/2019 By Michael Hackwill
%--------------------------------------------------------------------


%calculating array of sun vector positions at half hour incriments
%throughout the year
if isfield(meta,'Year')
    year=meta.Year;
else
    year=2013;
    meta.Year=2013;
end
sunvec=generate_sun_array(year,48,meta);

%calculating number of times to run through 
num=min([total,length(meta.solar_users)]);
if shadowonly
    num=min(num,length(shadowinglist));
end

%predefining outputs size
correlation=zeros(num,days,31);
res=zeros(num,4,2);
edges=zeros(num,days,31,2);
fprintf('house number ');

%begin loop for each house 
for i=1:num
    %finding the rolling correlation for the houshold
    if shadowonly
        itt=shadowinglist(i);
    else
        itt=i;
    end
    if mod(i,100)==99
        clc;
    end
    cap=s.solar_cap(itt);
    res(i,1,:)=[itt,itt];
    if noisered
        datai=zeros(days,48);
        for l=1:days
            for m=1:48
                temp=data(s.solar_users(itt),l,m);
                if (2*temp)>1
                    datai(l,m)=nthroot(2*temp,4)/2;
                else
                    datai(l,m)=temp;
                end
            end
        end
    else
        datai=squeeze(data(s.solar_users(itt),:,:));
    end
    fprintf('%3d\n',itt);
    %calculating sunrise/set timess 
    [sun,panel]=sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),0,s.dark_end,days);    
    %calculate 
    correlation=findcorrelation(i,itt,n,width,s,datai,meta,correlation,corrmanip,days);
    
    %predefining useful data
    correlate=squeeze(correlation(i,:,:));
    corrgrad=gradient(correlate);
    corrgrgrad=gradient(corrgrad);
    datai=squeeze(data(itt,:,10:40));
    datagradi=gradient(datai);
    
    %setting up edge detection for morning and afternoon; 
    datacap=-cap/8;
    morningshad=(((datai>datacap)&(correlate>correlationtrig))');
    afternoonshad=(((datai>datacap)&(correlate>correlationtrig))'); 
    %removing days that have no generation 
    for v=1:days
        if sum(squeeze(morningshad(:,v)))>28
            morningshad(:,v)=0;
            afternoonshad(:,v)=0;
        end 
    end
    morningshad(15:end,:)=0;            
    afternoonshad(1:15,:)=0;
    resmorn=qualitychecks(morningshad,showplot,sun,panel,s,1,days,meta);
    res(i,2:5,1)=resmorn;
    res(i,6,1)=datacap;
    res(i,7,1)=correlationtrig;
    resafter=qualitychecks(afternoonshad,showplot,sun,panel,s,2,days,meta);
    res(i,2:5,2)=resafter;
    res(i,6,2)=datacap;
    res(i,7,2)=correlationtrig;
    %clearing near sunrise and sunset
    width2=4;
    n2=20;
    for k=0:(days-1)
        for j=0:(31-width2)
            % w=data(s.solar_users(i),mod(((0:n)+k),days)+1,(1:3+j+s.dark_end));
            q=(0:(width2-1))+j+s.dark_end;
            temp1=sunvec((mod(((0:n2)+k),days)+1),q,3);
            a=double(squeeze(temp1));
            % v=s.solar_by_pc(s.postcode(i)==meta.pclist,mod(((0:n)+k),days)+1,1:3+j);
            temp2=correlate(mod(((0:n2)+k),days)+1,(1:width2)+j);
            b=squeeze(temp2);
            b(~any(~isnan(a), 2),:)=[];
            a(~any(~isnan(a), 2),:)=[];
            R=corr2(a,b);
            shift=floor(width2/2);
            correlation2(i,k+1,j+shift)=R;
        end
    end
    %quality checks
    %edge detection methods 
    [BW, thresholds]=edge(correlate','Canny');
    fprintf('%1.3f %1.3f\n',thresholds(1),thresholds(2));
    %removinig edges around sunrise 
    for w=1:31
        for q=1:days
            if (w<=(sun(1,q)-s.dark_end+2.5))||(w>=sun(2,q)-s.dark_end-2)
                BW(w,q)=0;
            end
            if (w<=(panel(1,q)-s.dark_end+2.5))||(w>=panel(2,q)-s.dark_end-2)
                BW(w,q)=0;
            end
        end
    end

    cc=bwconncomp(BW,8);
    shadstats=regionprops(cc);
    box=[shadstats.BoundingBox];
    box=reshape(box,[4,length(box)/4]);
    mornindex=find([shadstats.Area]>50&box(2,:)<12);
    afterindex=find([shadstats.Area]>50&(box(2,:)+box(4,:))>18);
    morningshadedge=ismember(labelmatrix(cc),mornindex);
    afternoonshadedge=ismember(labelmatrix(cc),afterindex);
    if sum(sum(morningshadedge))<200
        morningshadedge(:,:)=0;
    end
    if sum(sum(afternoonshadedge))<200
        afternoonshadedge(:,:)=0;
    end
    edges(i,:,:,1)=morningshadedge';
    edges(i,:,:,2)=afternoonshadedge';
    origafter=imdilate(afternoonshadedge,[0,0,0;0,1,0;1,1,1]);
    origmorn=imdilate(morningshadedge,[1,1,1;0,1,0;0,0,0]);
    for h=1:5
        origmorn=imdilate(origmorn,[1,1,1;0,1,0;0,0,0]);
        origafter=imdilate(origafter,[0,0,0;0,1,0;1,1,1]);        
    end
    lowmorndata=prctile(datai(origmorn'),5);
    highmorndata=prctile(datai(origmorn'),25);
    morn=datai((origmorn'&(datai>lowmorndata))&(datai<highmorndata));
    mornstd=log(std(morn)/cap);
    lowafterdata=prctile(datai(origafter'),5);
    highafterdata=prctile(datai(origafter'),25);
    after=datai((origafter'&(datai>lowafterdata))&(datai<highafterdata));
    afterstd=log(std(after)/cap);
    [~,mornerr]=find_sunvec(origmorn,morningshadedge,0,meta);
    [~,aftererr]=find_sunvec(origafter,afternoonshadedge,0,meta);
    zenith=sunvec(:,11:41,3)';
    azim=mod(sunvec(:,11:41,2)'+180,360);
    res(i,10:14,1)=[abs(1-mornerr),highmorndata,mornstd,thresholds];
    res(i,10:14,2)=[abs(1-aftererr),highafterdata,afterstd,thresholds];
    morningshadprob=mnrval(B,squeeze(res(i,[10,12,13],1)));
    afternoonshadprob=mnrval(B,squeeze(res(i,[10,12,13],2)));
    res(i,9,1)=morningshadprob(2);
    res(i,9,2)=afternoonshadprob(2);
    if showplot
        %plotting figures 
        figure(100);
        imagesc(datai(:,:)');
        hold on 
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end,days);        
        figure(101);
        imagesc(squeeze(correlation(i,:,:))')
        hold on
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end,days);
        axis([1 days 1 31])
        figure(102);
        imagesc(morningshadedge|afternoonshadedge);
        hold on
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end,days);
        axis([1 days 1 31])        
        figure(103);
        imagesc(squeeze(correlation2(i,:,:))');
        topdata=prctile(datai(origmorn'),10);
        fprintf('%3.5f\n',topdata);
        mornX=zenith(morningshadedge|afternoonshadedge);
        mornY=azim(morningshadedge|afternoonshadedge);
        figure(104);
        plot(mornY,mornX,'*');
        hold on
        plot(azim(1:end,170),zenith(1:end,170));
        hold on
        plot(azim(1:end,350),zenith(1:end,350));
        axis([60 300 0 100])
        set(gca, 'YDir','reverse');
        figure(105);
        imagesc(morningshadedge);
        figure(106);
        imagesc(afternoonshadedge);
        if vischeck
            res(i,9,1)=input('Morning shading True False\n');
            res(i,9,2)=input('Afternoon shading True False\n');
            fprintf(repmat('\b',1,72));
        end
        close([100 101 102 103 104 105 106]);
    end
fprintf(repmat('\b',1,16));
end
fprintf('\n')    
end

function correlation=findcorrelation(i,itt,n,width,s,datai,meta,correlation,corrmanip,days)
    [sun,panel]=sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),0,s.dark_end,days);
    for k=0:(days-1)
        for j=0:(31-width)
            % w=data(s.solar_users(i),mod(((0:n)+k),days)+1,(1:3+j+s.dark_end));
            q=(0:(width-1))+j+s.dark_end;
            temp1=datai(mod(((0:n)+k),days)+1,q);
            a=double(squeeze(temp1));
            % v=s.solar_by_pc(s.postcode(i)==meta.pclist,mod(((0:n)+k),days)+1,1:3+j);
            temp2=s.solar_by_pc(s.postcode(itt)==meta.pclist,mod(((0:n)+k),days)+1,(1:width)+j);
            b=squeeze(temp2);
            b(~any(~isnan(a), 2),:)=[];
            a(~any(~isnan(a), 2),:)=[];
            R=corr2(a,b);
            shift=floor(width/2);
            correlation(i,k+1,j+shift)=R;
        end
    end
    correlation=circshift(correlation,ceil(n/2),2);
    if corrmanip
        for m=0:364
            for n=0:(31-width)
                R=correlation(i,m+1,n+shift);
                if (n<=(panel(1,m+1)-s.dark_end-shift))||(n>=panel(2,m+1)-s.dark_end-shift)
                    correlation(i,m+1,n+shift)=abs(R+1)^(1/2)-1;%weaken correlation
                end
                if (n<=(sun(1,m+1)-s.dark_end-shift+1))||(n>=sun(2,m+1)-s.dark_end-shift)
                        correlation(i,m+1,n+shift)=0; %(1,2);
                else
                        correlation(i,m+1,n+shift)=R;
                end
            end   
        end
    end
end

function res=qualitychecks(shading,showplot,sun,panel,s,ampm,days,meta)
    trac=0;
    orig=shading;
    shading1=shading;
    %fix the shading 
    for w=1:31
        for q=1:days
            if (w<=(sun(1,q)-s.dark_end+1))||(w>=sun(2,q)-s.dark_end-1)
                shading(w,q)=0;
                shading1(w,q)=1;
            end
            if (w<=(panel(1,q)-s.dark_end+1))||(w>=panel(2,q)-s.dark_end-1)
                shading(w,q)=0;
                shading1(w,q)=1;
            end
        end
    end
    if mod(ampm,2)
        shading(15:end,:)=0;
        shading1(15:end,:)=0;
    else
        shading(1:15,:)=0;
        shading1(1:15,:)=0;
    end    
    for q=1:days
        if sum(shading(:,q))>=1
            trac=trac+1;
        end
    end
    

        
    
    %first condition
    res(1)=trac/days;
    
    %condition 2
    shading=shading1;
    shading=bwareaopen(shading,5);
    shading=imdilate(shading,[0 1 0;1 1 1;0 1 0]);
    shading=imabsdiff(shading, shading1);
    res(2)=360/length(find(shading));
    % third contion needs optimisation and fixing
    if res(1)>0.5
    shading=orig;
    for j=1:2
        shading=imerode(shading,[1 1 1;1 1 1;1 1 1]);
    end
    shading=imabsdiff(shading, orig);
    for w=1:31
        for q=1:days
            if (w<=(sun(1,q)-s.dark_end+1))||(w>=sun(2,q)-s.dark_end-1)
                shading(w,q)=0;
            end
            if (w<=(panel(1,q)-s.dark_end+1))||(w>=panel(2,q)-s.dark_end-1)
                shading(w,q)=0;
            end
        end
    end
    [~,err]= find_sunvec(orig,shading,showplot,meta);
    res(3)=err;
    else
        res(3)=0;
    end
    
    res(4)=(res(1)>0.8)&(res(2)>0.65)&(res(3)>0.7);
end


    
    
    

    