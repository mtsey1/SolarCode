function [correlation,res]=cloudcorr(s,data,meta) 
n=60; %length of correlation
width=3; % width of correlation
total=100;
showplot=0;
shadowonly=0;
corrmanip=1;
gencheck=0;
noisered=1;
modelfit=0;
shadowinglist=[20,146,179,208,233,238,239,240,243,259,266,277,284,292];
%shadowinglist=[15 19 20 22 25 27 31 35 37 42 46 54 56 58 60 64 69 71 76 77 84 88 89];
%shadowinglist=[18 25 31 64 76 89 98];

modelfun=@(c,x) c(1)+c(2)*cosd(360/365*x+c(3)); 

num=total;
if shadowonly
    num=min(total,length(shadowinglist));
end


correlation=zeros(num,365,31);
shading=zeros(2,num);
mornshadingcoeff=zeros(num,3);
aftershadingcoeff=zeros(num,3);
res=zeros(num,4,2);
for i=1:num
    %finding the rolling correlation for the houshold
    if shadowonly
        itt=shadowinglist(i);
    else
        itt=i;
    end
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
        imagesc(datai');
        figure(2);
    end
    imagesc(squeeze(data(s.solar_users(itt),:,:))');
    [sun,panel]=sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),0,s.dark_end);    
    correlation=findcorrelation(i,itt,n,width,s,datai,meta,correlation,corrmanip);
    %{
    if shadowonly
        itt=shadowinglist(i);
    else
        itt=i;
    end
    [sun,panel]=sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),0,s.dark_end);
    for k=0:364
        for j=0:(31-width)
            % w=data(s.solar_users(i),mod(((0:n)+k),365)+1,(1:3+j+s.dark_end));
            q=(0:(width-1))+j+s.dark_end;
            temp1=data(s.solar_users(itt),mod(((0:n)+k),365)+1,q);
            a=double(squeeze(temp1));
            % v=s.solar_by_pc(s.postcode(i)==meta.pclist,mod(((0:n)+k),365)+1,1:3+j);
            temp2=s.solar_by_pc(s.postcode(itt)==meta.pclist,mod(((0:n)+k),365)+1,(1:width)+j);
            b=squeeze(temp2);
            %R=corr(a',b');
            R=corr2(a,b);
            shift=floor(width/2);
            if (j<=(panel(1,k+1)-s.dark_end-shift))||(j>=panel(2,k+1)-s.dark_end-shift)
                correlation(i,k+1,j+shift)=abs(R)^(1/2)*R;%weaken correlation
            end
            if (j<=(sun(1,k+1)-s.dark_end-shift))||(j>=sun(2,k+1)-s.dark_end-shift)
                    correlation(i,k+1,j+shift)=0; %(1,2);
            else
                    correlation(i,k+1,j+shift)=R;
            end
                
        end
    end
    %}
    %predefining useful data
    cap=s.solar_cap(itt);
    post=s.postcode(itt);
    correlate=squeeze(correlation(i,:,:));
    corrgrad=gradient(correlate);
    corrgrgrad=gradient(corrgrad);
    datai=squeeze(data(itt,:,10:40));
    datagradi=gradient(datai);
    
    %setting up edge detection for morning and afternoon; 
    morningshad=(((datai>-cap/5)&(correlate>-0.25))');
    morningshad(15:end,:)=0;
    afternoonshad=(((datai>-cap/5)&(correlate>-0.25))'); 
    afternoonshad(1:15,:)=0;
    
    %clearing near sunrise and sunset
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
    
    fprintf('house number %d\n',itt);
    
    %quality checks 
    res=qualitychecks(i,res,morningshad,1,showplot);
    res=qualitychecks(i,res,afternoonshad,2,showplot);
    
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
    if showplot
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
        figure(101); 
        %imagesc((datai>0)')
        %imagesc(imgaussfilt(squeeze(data(i,:,10:41))',[0.5 2]));
        imagesc(squeeze(data(itt,:,10:41))');
        hold on
        if shading(1,i)
            plot(modelfun(mornshadingcoeff(i,:),1:365))
            hold on
        end
        if shading(2,i)
            plot(modelfun(aftershadingcoeff(i,:),1:365))
            hold on
        end
        figure(102);
        imagesc(afternoonshad)
        hold on
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end);
        axis([1 365 1 31])
        %imagesc(squeeze(s.solar_by_pc(post==meta.pclist,:,:))');
        %{
        figure(1);
        plot(squeeze(data(i,10,:)));
        figure(2);
        plot(squeeze(data(i,178,:)));
        %}
        figure(103);
        imagesc(corrgrad')
        figure(105);
        imagesc(corrgrgrad');
        figure(104);
        imagesc(morningshad);    %hold on;
        %sunriseset(meta.location.latitude,meta.location.longitude,11,s.solar_az(i),(s.solar_ze(i)),1,s.dark_end);
        %imagesc((squeeze(correlation(i,:,:))'>-0.2))
        hold on
        sunriseset(meta.location.latitude,meta.location.longitude,10,s.solar_az(itt),(s.solar_ze(itt)),1,s.dark_end);
        axis([1 365 1 31])
        %count(i)=nnz(squeeze((correlation(i,:,:)>-0.2)&(correlation(i,:,:)~=0))')
        %}
        close([100 101 104 102]);
    end
end
        
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
            R=corr2(a,b);
            shift=floor(width/2);
            correlation(i,k+1,j+shift)=R;
            if corrmanip
                if (j<=(panel(1,k+1)-s.dark_end-shift))||(j>=panel(2,k+1)-s.dark_end-shift)
                    correlation(i,k+1,j+shift)=abs(R)^(1/2)*R;%weaken correlation
                end
                if (j<=(sun(1,k+1)-s.dark_end-shift))||(j>=sun(2,k+1)-s.dark_end-shift)
                        correlation(i,k+1,j+shift)=0; %(1,2);
                else
                        correlation(i,k+1,j+shift)=R;
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

function res=qualitychecks(i,res,shading,ampm,showplot)
    trac=0;
    for q=1:365
        if sum(shading(:,q))>=1
            trac=trac+1;
        end
    end
    res(i,2,ampm)=trac/365;
    perim=length(find(bwperim(shading),8));
    res(i,3,ampm)=2*trac/perim;
    [~,err]=find_sunvec(shading,showplot);
    res(i,4,ampm)=err;
end
