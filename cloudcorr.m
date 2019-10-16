function [correlation,res,edges,meta,fitvals,simed]=cloudcorr(s,data,meta) 
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
total=300; %Max number of household analysed  
showplot=0; % logical value as to whether to show plots and figures generated 
shadowonly=0; %logical to use subset of data specified by vector shadowinglist
corrmanip=1; %logical to manipulate correlation plot around sunrise/set
noisered=1; %logical to reduce the amplitude of usage noise
vischeck=1; %logical value allowing to input the
res=1000;
%shadowinglist: a list of specific houses that are needed to be checked 
%shadowinglist=[13,14,15,19,24,30,146,179,208,233,238,239,240,243,259,266,277,284,292];
shadowinglist=[5,21,37,53,69,85,101];
%UNUSED-------------------------------------------------------
%B: old logistic regrestion values
B=[10.8481630139160;-16.2828877554926;-1.06561285029183;5.08033613741109];
correlationtrig=-0.35; %threshold //Curretly unused
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
%Edited 16/10/2019 By Michael Hackwill
%--------------------------------------------------------------------


%calculating array of sun vector positions at half hour incriments
%throughout the year

%filling year field
if isfield(meta,'Year')
    year=meta.Year;
else
    year=2013;
    meta.Year=2013;
end


%creating arrays for sun position ove the course of the year
sunvec=generate_sun_array(year,48,meta);
sunarray=generate_sun_array(year,res,meta);

%calculating number of times to run through 
num=min([total,length(meta.solar_users)]);
if shadowonly
    num=min(num,length(shadowinglist));
end

%predefining outputs size
correlation=zeros(num,days,31);
res=zeros(num,4,2);
edges=zeros(num,days,31,2);
fitvals=zeros(num,2,2);
simed=[];
fprintf('house number ');

%begin loop for each house 
for i=1:num
    
    %extracting solar users 
    if shadowonly
        itt=shadowinglist(i);
    else
        itt=i;
    end
    
    %clearing output every 100 houses 
    if mod(i,100)==99
        clc;
    end
    
    %extracting elsimated solar capacity of the household
    cap=s.solar_cap(itt);
    res(i,1,:)=[itt,itt];
    
    
    
    if noisered
        %-----------------------------------------
        %This section attempts to reduce the effect of consumption noise by
        %trimming the peak value by taking the 4th root of each value above
        %0.5
        % I admit this is an arbitary work around and a better method may
        % yeild better results although this does improve the outputs
        % slightly
        %-------------------------------------------------------------
        datai=zeros(days,48);
        %checking over full array of generation 
        for l=1:days
            for m=1:48
                temp=data(s.solar_users(itt),l,m);
                if (2*temp)>1
                    datai(l,m)=nthroot(temp*4,4)/4;
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
    %calculate correlation array from rolling correlation function 
    correlation=findcorrelation(i,itt,n,width,s,datai,meta,correlation,corrmanip,days);
    
    %predefining useful data
    correlate=squeeze(correlation(i,:,:));
    datai=squeeze(datai(itt,:,10:40));
    %UNUSED-----------------
    corrgrad=gradient(correlate);
    corrgrgrad=gradient(corrgrad);
    datagradi=gradient(datai);
    %--------------------------



    
    %----------------------------------------------------
    %-----------------Edge Detection--------------------
    %---------------------------------------------------
    %Edge detection is one region which can be improved current method uses
    %basic gradient methods to find edges. This is one area where a more
    %robust method could be developed/implemented for vastly superior results as  
    %setting up edge detection for morning and afternoon; 
    
    %datacap sets value for which an edge may be discarded for being above
    %to red
    datacap=-cap/8;
    
    %UNUSED------------------------------------------------------
    morningshad=(((datai>datacap)&(correlate>correlationtrig))');
    afternoonshad=(((datai>datacap)&(correlate>correlationtrig))');
    %---------------------------------------------------------------

    
    
    %++++++++++++++++++++++++++++OLD EDGE METHOD++++++++++++++++++++
    %removing days that have no generation 
    %attempted work around for days when solar panel is switched off
    %probably redunant so could be removed for optimisation
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
    %corelation between correlation and the solar angle to try to see edge
    %again can be probably got rid of
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
    %++++++++++++++++++++++++++++END OF OLD METHOD+++++++++++++++++++++
    %edge detection methods using Canny see MATLAB guide for justification  
    [BW, thresholds]=edge(correlate','Canny');
    %fprintf('%1.3f %1.3f\n',thresholds(1),thresholds(2));
    %removinig any edges that may possibly be caused by sunrise rather than
    %an object
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

    %checking for connected edges
    cc=bwconncomp(BW,8);
    
    %extracting states of the connected regions 
    shadstats=regionprops(cc);
    
    %find indexes of the box that bound the edges to get rid of small edges
    %and seperate the edges into morning and afternoon edges 
    box=[shadstats.BoundingBox];
    box=reshape(box,[4,length(box)/4]);
    
    %seperating edges 
    mornindex=find([shadstats.Area]>50&box(2,:)<12);
    afterindex=find([shadstats.Area]>50&(box(2,:)+box(4,:))>18);
    morningshadedge=ismember(labelmatrix(cc),mornindex);
    afternoonshadedge=ismember(labelmatrix(cc),afterindex);
    
    %ignoring edges when they dont consist at least along 200 days
    if sum(sum(morningshadedge))<200
        morningshadedge(:,:)=0;
    end
    if sum(sum(afternoonshadedge))<200
        afternoonshadedge(:,:)=0;
    end
    
    %enteringdata into array 
    edges(i,:,:,1)=morningshadedge';
    edges(i,:,:,2)=afternoonshadedge';
    
    %shifting edges back to counteract both the shift due to correlation
    %and edge detection errors
    edges(i,:,:,1)=circshift(squeeze(edges(i,:,:,1)),-2,2);
    edges(i,:,:,2)=circshift(squeeze(edges(i,:,:,2)),2,2);
    gen=datai>0;
    edges(i,:,:,1)=(squeeze(edges(i,:,:,1))&gen);
    edges(i,:,:,2)=(squeeze(edges(i,:,:,2))&gen);
    %could potentially remove any tiny points from
    
    
    %below can be replaced with post fit line data check
    origafter=imdilate(afternoonshadedge,[0,0,0;0,1,0;1,1,1]);
    origmorn=imdilate(morningshadedge,[1,1,1;0,1,0;0,0,0]);
    for h=1:5
        origmorn=imdilate(origmorn,[1,1,1;0,1,0;0,0,0]);
        origafter=imdilate(origafter,[0,0,0;0,1,0;1,1,1]);        
    end
    
    
    %check the generation in the shadded region to see if there is a large amount of generation in the shaded region
    
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
    [simed,fitvals]=simedge(edges,meta,i,sunvec,sunarray,showplot,fitvals,simed);
    
    
    if showplot
        %plotting figures 
        mornedg=squeeze(edges(i,:,:,1));
        afteredg=squeeze(edges(i,:,:,2));
        
        %Plot of trimmed data with sunrise and sunset 
        figure(100);
        imagesc(datai(:,:)');
        hold on 
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end,days);        
        
        %plot of correlation with sunrise and sunset
        figure(101);
        correlation(i,squeeze(correlation(i,:,:)>0))=0;
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
        row=simed(i,1,:);
        col=simed(i,2,:);
        %figure(107);
        %imagesc(squeeze(data(i,:,10:40))'); hold on;
        %plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
        %figure(108);
        %imagesc(squeeze(correlation(i,:,:))'); hold on;
        %plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
        %figure(109);
        %imagesc((mornedg|afteredg)');hold on;
        %plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
        %axis([1 days 1 31])
        if vischeck
            res(i,9,1)=input('Morning shading True False\n');
            res(i,9,2)=input('Afternoon shading True False\n');
            fprintf(repmat('\b',1,72));
        end
        close([100 101 102 103 104 105 106 110]);
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
    correlation(i,:,:)=circshift(correlation(i,:,:),ceil(n/2),2);
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

function [simed,fitvals]=simedge(edges,meta,k,sunvec,sunarray,showplot,fitvals,simed)
year=meta.Year;
res=1000;
days=365+~mod(year,4);
mesh=zeros(days, res);
%sunarray=generate_sun_array(year,res,meta);
sunarray(:,:,2)=mod(sunarray(:,:,2)+180,360);
%{
leftheight=66;
rightheight=66;
leftedge=110;
rightedge=190;
centerheight=66;
%}
morn=0;
after=0;
sunvec2=sunvec;
%sunvec2(:,:,2)=mod(sunvec2(:,:,2)+180,360);
mornedg=squeeze(edges(k,:,:,1));
afteredg=squeeze(edges(k,:,:,2));
zenith2=squeeze(sunvec2(:,11:41,3)');
azim2=squeeze(mod(sunvec2(:,11:41,2)'+180,360));
mornzen=zenith2(logical(mornedg'));
mornazim=azim2(logical(mornedg'));
afterzen=zenith2(logical(afteredg'));
afterazim=azim2(logical(afteredg'));
if length(find(mornedg))>100
mornfit=fit(mornazim,mornzen,'poly1','Robust','on'); 
morn=1;
end
if length(find(afteredg))>100
afterfit=fit(afterazim,afterzen,'poly1','Robust','on'); 
after=1;
end
%afterfit=polyfit(afterazim,afterzen,1);

%----
%removingoutliers




for i=1:days
    for j=1:res
        if sunarray(i,j,3)>90 
            mesh(i,j)=1;
        %elseif polyval(pmorn,sunarray(i,j,2))<=sunarray(i,j,3)&&pmorn(1)<0&&sunarray(i,j,2)>180
        elseif morn
            if mornfit(sunarray(i,j,2))<=sunarray(i,j,3)&&sunarray(i,j,2)>180
            mesh(i,j)=1;
            end
        %elseif polyval(pafter,sunarray(i,j,2))<=sunarray(i,j,3)&&pafter(1)>0&&sunarray(i,j,2)<180
        end
        if after
            if afterfit(sunarray(i,j,2))<=sunarray(i,j,3)&&sunarray(i,j,2)<180
            mesh(i,j)=1;
            end
        %if sunarray(i,j,3)>90 
            %mesh(i,j)=1;
        %elseif sunarray(i,j,2)<=leftedge && sunarray(i,j,3)>leftheight
            %mesh(i,j)=1;
        %elseif sunarray(i,j,2)>=rightedge && sunarray(i,j,3)>rightheight
            %mesh(i,j)=1;
        %elseif (sunarray(i,j,2)<rightedge && sunarray(i,j,2)>leftedge) && sunarray(i,j,3)>centerheight 
            %mesh(i,j)=1;
        end 
    end
end




edge=diff(mesh');



%look to seperate morning and afternoon edge
[col,row]=find(edge);
col=col/1000*48-9;
simed(k,1,1:length(row))=row';
simed(k,2,1:length(col))=col';
if morn
    fitvals(k,:,1)=coeffvalues(mornfit);
else
    fitvals(k,:,1)=[0,90];
end
if after
    fitvals(k,:,2)=coeffvalues(afterfit);
else
    fitvals(k,:,2)=[0,90];
end
%meshmorn=zeros(365,31);
%meshafter=zeros(365,31);
%{
for p=1:days
    for l=10:40
        %elseif polyval(pmorn,sunarray(i,j,2))<=sunarray(i,j,3)&&pmorn(1)<0&&sunarray(i,j,2)>180
        if morn
            if mornfit(sunvec2(p,l,2))<=sunvec2(p,l,3)&&sunvec2(p,l,2)>180
            meshmorn(p,l-9)=1;
            end
        %elseif polyval(pafter,sunarray(i,j,2))<=sunarray(i,j,3)&&pafter(1)>0&&sunarray(i,j,2)<180
        end
        if after
            if afterfit(sunvec2(p,l,2))<=sunvec2(p,l,3)&&sunvec2(p,l,2)<180
            meshafter(p,l-9)=1;
            end
        end
    end
end
%}



%x=60:300;
%Ymorn=polyval(pmorn,x);
%Yafter=polyval(pafter,x);
%X=[60, leftedge, leftedge, rightedge, rightedge,300];
%Y=[leftheight,leftheight,centerheight,centerheight,rightheight,rightheight];

if showplot
zenith=sunarray(:,150:800,3)';
azim=sunarray(:,150:800,2)';
simed(k,:,:)=[row; col];
mornX=zenith2((mornedg|afteredg)');
mornY=azim2((mornedg|afteredg)');
figure(110);
plot(mornY,mornX,'*');
hold on
plot(azim(1:end,170),zenith(1:end,170));
hold on
plot(azim(1:end,350),zenith(1:end,350));
hold on
if morn
plot(mornfit);
hold on
end
if after
plot(afterfit);
end
axis([60 300 0 100])
set(gca, 'YDir','reverse');

%figure
%imagesc(mesh');
%current problem is when only one sided 
% problem with ouliers
%set(gca, 'YDir','reverse')
end

end

    
    
    

    