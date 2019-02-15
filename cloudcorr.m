function [correlation,res]=cloudcorr(s,data,meta) 
n=60; %length of correlation
width=3; % width of correlation
numb=1;
correlationtrig=linspace(-0.35,-0.15,numb);
days=365;
%correlationtrig=-0.3;
total=300;
showplot=1;
shadowonly=1;
corrmanip=1;
gencheck=0;
noisered=1;
modelfit=0;
removeriseset=0;
qualityparam=1;
vischeck=1;
B=[10.8481630139160;-16.2828877554926;-1.06561285029183;5.08033613741109];
shadowinglist=[20,146,179,208,233,238,239,240,243,259,266,277,284,292];
%shadowinglist=[157,220,23,213,182,244,299,264,243,204,80,208,26,152,20,134,16,94,55,114,119,4];
%shadowinglist=[291,197,40,108,103,84,173,20,213,277,56];
%shadowinglist=[15 19 20 22 25 27 31 35 37 42 46 54 56 58 60 64 69 71 76 77 84 88 89];
%shadowinglist=[18 25 31 64 76 89 98];
%shadowinglist=[1 5 7 8 16 18 20 22 29 30 33 35 39 41 42 43 44 45 46 50];
%afternoon shading
%shadowinglist=[89 92 64 39 70 33 8 95 34 17 19 43 72 20 45 52 31 46 1 6 58 75 91 23 51 17 26];
%mornign shading 
%shadowinglist=[24 31 55 88 8 100 35 20 81 94 52 25 79 30 44 41 69 16 60 74 90 42 98 23 80 58 93 18 57 56];
%shadowinglist=[241,233,287,154,252,42,265,269,236,212,78,295,140];
modelfun=@(c,x) c(1)+c(2)*cosd(360/365*x+c(3)); 
sunvec=generate_sun_array(2013,48);

num=total;
if shadowonly
    num=min(total,length(shadowinglist));
end


correlation=zeros(num,365,31);
shading=zeros(2,num);
mornshadingcoeff=zeros(num,3);
aftershadingcoeff=zeros(num,3);
res=zeros(num,4,2);
fprintf('house number ');   
for i=1:num
    %finding the rolling correlation for the houshold
    if shadowonly
        itt=shadowinglist(i);
    else
        itt=i;
    end
    cap=s.solar_cap(itt);
    res(i,1,:)=[itt,itt];
    if noisered
        datai=zeros(365,48);
        for l=1:365
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
    if showplot
        figure(1);
        imagesc(datai(:,10:41)');
        figure(2);
        imagesc(squeeze(data(s.solar_users(itt),:,10:40))');
    end
    fprintf('%3d\n',itt);
    [sun,panel]=sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),0,s.dark_end);    
    correlation=findcorrelation(i,itt,n,width,s,datai,meta,correlation,corrmanip);
    
    %predefining useful data
    post=s.postcode(itt);
    correlate=squeeze(correlation(i,:,:));
    corrgrad=gradient(correlate);
    corrgrgrad=gradient(corrgrad);
    datai=squeeze(data(itt,:,10:40));
    datagradi=gradient(datai);
    
    %setting up edge detection for morning and afternoon; 
    thrsh=cap/8;
    datacap=linspace(-cap/8-thrsh,-cap/8+thrsh,numb);
    mornmax=0;
    aftermax=0;
    for j=1:numb
        for k=1:numb
            morningshad=(((datai>datacap(j))&(correlate>correlationtrig(k)))');
            afternoonshad=(((datai>datacap(j))&(correlate>correlationtrig(k)))'); 
            %probably more efficient way of doing this 
            for days=1:365
                if sum(squeeze(morningshad(:,days)))>28
                    morningshad(:,days)=0;
                    afternoonshad(:,days)=0;
                end 
            end
            morningshad(15:end,:)=0;            
            afternoonshad(1:15,:)=0;
            if qualityparam
                resmorn=qualitychecks(morningshad,showplot,sun,panel,s,1);
                %may need to be adjusted 
                summorn=sum(resmorn(1:3));
                if summorn>mornmax
                    mornmax=summorn;
                    res(i,2:5,1)=resmorn;
                    res(i,6:7,1)=[datacap(j),j];
                    res(i,8:9,1)=[correlationtrig(k),k];
                end
                resafter=qualitychecks(afternoonshad,showplot,sun,panel,s,2);
                sumafter=sum(resafter(1:3));
                if sumafter>aftermax
                    aftermax=sumafter;
                    res(i,2:5,2)=resafter;
                    res(i,6:7,2)=[datacap(j),j];
                    res(i,8:9,2)=[correlationtrig(k),k];
                end
            end
        end
    end
    %clearing near sunrise and sunset
    if removeriseset
        for w=1:31
            for q=1:365
            if (w<=(sun(1,q)-s.dark_end+1))||(w>=sun(2,q)-s.dark_end-1)
                morningshad(w,q)=0;
                afternoonshad(w,q)=0;
            end
            if (w<=(panel(1,q)-s.dark_end+1))||(w>=panel(2,q)-s.dark_end-1)
                    morningshad(w,q)=0;
                    afternoonshad(w,q)=0;
            end
        end
        end
    end
    width2=4;
    n2=20;
    for k=0:364
        for j=0:(31-width2)
            % w=data(s.solar_users(i),mod(((0:n)+k),365)+1,(1:3+j+s.dark_end));
            q=(0:(width2-1))+j+s.dark_end;
            temp1=sunvec((mod(((0:n2)+k),365)+1),q,3);
            a=double(squeeze(temp1));
            % v=s.solar_by_pc(s.postcode(i)==meta.pclist,mod(((0:n)+k),365)+1,1:3+j);
            temp2=correlate(mod(((0:n2)+k),365)+1,(1:width2)+j);
            b=squeeze(temp2);
            b(~any(~isnan(a), 2),:)=[];
            a(~any(~isnan(a), 2),:)=[];
            R=corr2(a,b);
            shift=floor(width2/2);
            correlation2(i,k+1,j+shift)=R;
        end
    end
    %quality checks 

    % need to add delay to the model as currently miss aligned?
    %morning shading
    if modelfit
        [shading,mornshadingcoeff,aftershadingcoeff]=modelfitt(shading,i,morningshad,1,200,0.5,mornshadingcoeff,aftershadingcoeff);
        %{
        ii=find(morningshad);
        [yy,xx]=ind2sub(size(morningshad),ii);
        if length(unique(xx))<200
            %maybe can change from size to to the number of unique days
            fprintf('no morning shading\n')
            shading(1,i)=0;
        else
            mdl=fitnlm(xx,yy,modelfun,[9,3,-81]);
            if mdl.Rsquared.adjusted<0.5
                % maybe rmse would be a better gauge 
                shading(1,i)=0;
                fprintf('no morning shading\n')
            else
                shading(1,i)=1;
                fprintf('morning shading\n')
                mornshadingcoeff(i,:)=mdl.Coefficients{:,'Estimate'};
            end
        end 
        %}
        %afternoon shading 
        [shading,mornshadingcoeff,aftershadingcoeff]=modelfitt(shading,i,afternoonshad,2,200,0.5,mornshadingcoeff,aftershadingcoeff);    
        %{
        ii=find(afternoonshad);
        [yy,xx]=ind2sub(size(afternoonshad),ii);
        if length(unique(xx))<200
            fprintf('no afternoon shading\n')
            shading(2,i)=0;
        else    
            mdl=fitnlm(xx,yy,modelfun,[25,3,-81]);
            if mdl.Rsquared.adjusted<0.5
                shading(2,i)=0;
                fprintf('no afternoon shading \n')
            else
                shading(2,i)=1;
                fprintf('afternoon shading \n')
                aftershadingcoeff(i,:)=mdl.Coefficients{:,'Estimate'};
            end
        end
        %}
    end
    if gencheck
        if shading(1,i)
            for r=1:365
                v=floor(modelfun(mornshadingcoeff(i,:),r));
                for p=1:v
                    if (datai(r,p)<-cap/5)
                        shading(1,i)=0;
                        break;
                    end
                end
                if ~shading(1,i)
                    fprintf('no morning shading\n')
                    break;
                end
            end
        end
        if shading(2,i)
            for k=1:365
                m=ceil(modelfun(aftershadingcoeff(i,:),k));
                for h=m:30
                    if (datai(k,h+1)<-cap/5)
                        shading(2,i)=0;
                        break;
                    end
                end
                if ~shading(2,i)
                    fprintf('no afternoon shading\n')
                    break;
                end
            end
        end
    end
    %Maybe test eliminating outliers before fitting or after and iterate fit again 
    %figure plotting 
    [BW, thresholds]=edge(correlate','Canny');
    fprintf('%1.3f %1.3f\n',thresholds(1),thresholds(2));
    for w=1:31
        for q=1:365
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
    [~,mornerr]=find_sunvec(origmorn,morningshadedge,0);
    [~,aftererr]=find_sunvec(origafter,afternoonshadedge,0);
    zenith=sunvec(:,11:41,3)';
    azim=mod(sunvec(:,11:41,2)'+180,360);
    res(i,10:14,1)=[abs(1-mornerr),highmorndata,mornstd,thresholds];
    res(i,10:14,2)=[abs(1-aftererr),highafterdata,afterstd,thresholds];
    morningshadprob=mnrval(B,squeeze(res(i,[10,12,13],1)));
    afternoonshadprob=mnrval(B,squeeze(res(i,[10,12,13],2)));
    res(i,9,1)=morningshadprob(2);
    res(i,9,2)=afternoonshadprob(2);
    if showplot
        morningshad=(((datai>res(i,6,1))&(correlate>res(i,8,1)))');
        afternoonshad=(((datai>res(i,6,2))&(correlate>res(i,8,2)))'); 
        figure(100);
        imagesc(squeeze(correlation(i,:,:))')
        hold on
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end);
        hold on 
        if shading(1,i)
            plot(modelfun(mornshadingcoeff(i,:),1:365))
            hold on
        end
        if shading(2,i)
            plot(modelfun(aftershadingcoeff(i,:),1:365))
            hold on
        end
        axis([1 365 1 31])
        %figure(101); 
        %imagesc((datai>0)')
        %imagesc(imgaussfilt(squeeze(data(i,:,10:41))',[0.5 2]));
        %{
        hold on
        if shading(1,i)
            plot(modelfun(mornshadingcoeff(i,:),1:365))
            hold on
        end
        if shading(2,i)
            plot(modelfun(aftershadingcoeff(i,:),1:365))
            hold on
        end
        ---------------------------------------------
        figure(102);
        imagesc(afternoonshad)
        hold on
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end);
        axis([1 365 1 31])
        %imagesc(squeeze(s.solar_by_pc(post==meta.pclist,:,:))');
        ---------------------------------------------
        figure(1);
        plot(squeeze(data(i,10,:)));
        figure(2);
        plot(squeeze(data(i,178,:)));
        %}
        if 1==0
            figure(103);
            imagesc(corrgrad')
            figure(105);
            imagesc(corrgrgrad');
        
        figure(104);
        imagesc(morningshadedge|afternoonshadedge);
        %hold on;
        %sunriseset(meta.location.latitude,meta.location.longitude,11,s.solar_az(i),(s.solar_ze(i)),1,s.dark_end);
        %imagesc((squeeze(correlation(i,:,:))'>-0.2))
        hold on
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end);
        axis([1 365 1 31])
        end
        %count(i)=nnz(squeeze((correlation(i,:,:)>-0.2)&(correlation(i,:,:)~=0))')
        %}
        
        figure(106);
        imagesc(squeeze(correlation2(i,:,:))');
        topdata=prctile(datai(origmorn'),10);
        fprintf('%3.5f\n',topdata);
        mornX=zenith(morningshadedge|afternoonshadedge);
        mornY=azim(morningshadedge|afternoonshadedge);
        figure(107);
        plot(mornY,mornX,'*');
        hold on
        plot(azim(1:end,170),zenith(1:end,170));
        hold on
        plot(azim(1:end,350),zenith(1:end,350));
        axis([60 300 0 100])
        set(gca, 'YDir','reverse');

        %hold on;
        %plot(x,y);
        %{
        figure(108);
        afterX=zenith(afternoonshadedge);
        afterY=azim(afternoonshadedge);
        plot(afterY,afterX,'*')
        %}
        figure(109);
        imagesc(morningshadedge);
        figure(110);
        imagesc(afternoonshadedge);
        if vischeck
            res(i,9,1)=input('Morning shading True False\n');
            res(i,9,2)=input('Afternoon shading True False\n');
            fprintf(repmat('\b',1,72));
        end
        close([100 101 104 102 107 108]);
    end
fprintf(repmat('\b',1,16));
end
fprintf('\n')    
end

function correlation=findcorrelation(i,itt,n,width,s,datai,meta,correlation,corrmanip)
    [sun,panel]=sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),0,s.dark_end);
    for k=0:364
        for j=0:(31-width)
            % w=data(s.solar_users(i),mod(((0:n)+k),365)+1,(1:3+j+s.dark_end));
            q=(0:(width-1))+j+s.dark_end;
            temp1=datai(mod(((0:n)+k),365)+1,q);
            a=double(squeeze(temp1));
            % v=s.solar_by_pc(s.postcode(i)==meta.pclist,mod(((0:n)+k),365)+1,1:3+j);
            temp2=s.solar_by_pc(s.postcode(itt)==meta.pclist,mod(((0:n)+k),365)+1,(1:width)+j);
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

function [shading,mornshadcoeff,aftershadcoeff]=modelfitt(shading,i,shadeedge,ampm,daysthresh,corrthresh,mornshadcoeff,aftershadcoeff)
modelfun=@(c,x) c(1)+c(2)*cosd(360/365*x+c(3));    
ii=find(shadeedge);
    [yy,xx]=ind2sub(size(shadeedge),ii);
    if length(unique(xx))<daysthresh
        %maybe can change from size to to the number of unique days
            if mod(ampm,2)
                fprintf('no morning shading\n')
            else
                fprintf('no afternoon shading\n')
            end
        shading(ampm,i)=0;
    else
        mdl=fitnlm(xx,yy,modelfun,[9,3,-81]);
        if mdl.Rsquared.adjusted<corrthresh
            % maybe rmse would be a better gauge 
            shading(ampm,i)=0;
            if mod(ampm,2)
                fprintf('no morning shading\n')
            else
                fprintf('no afternoon shading\n')
            end
        else
            shading(ampm,i)=1;
            if mod(ampm,2)
                fprintf('morning shading\n')
                mornshadcoeff(i,:)=mdl.Coefficients{:,'Estimate'};
            else
                fprintf('afternoon shading\n')
                aftershadcoeff(i,:)=mdl.Coefficients{:,'Estimate'};
            end
        end
    end
end

function res=qualitychecks(shading,showplot,sun,panel,s,ampm)
    trac=0;
    orig=shading;
    shading1=shading;
    %fix the shading 
    for w=1:31
        for q=1:365
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
    for q=1:365
        if sum(shading(:,q))>=1
            trac=trac+1;
        end
    end
    

        
    
    %first condition
    res(1)=trac/365;
    
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
        for q=1:365
            if (w<=(sun(1,q)-s.dark_end+1))||(w>=sun(2,q)-s.dark_end-1)
                shading(w,q)=0;
            end
            if (w<=(panel(1,q)-s.dark_end+1))||(w>=panel(2,q)-s.dark_end-1)
                shading(w,q)=0;
            end
        end
    end
    [~,err]=find_sunvec(orig,shading,showplot);
    res(3)=err;
    else
        res(3)=0;
end
    
    res(4)=(res(1)>0.8)&(res(2)>0.65)&(res(3)>0.7);
end


    
    
    

    