function [simed,fitvals]=simedge(corr,edges,data,meta,p)
close all;
year=meta.Year;
res=1000;
days=365+~mod(year,4);
mesh=zeros(days, res);
sunarray=generate_sun_array(year,res,meta);
sunarray(:,:,2)=mod(sunarray(:,:,2)+180,360);
morn=0;
after=0;
all=1;
preload=1;
sunvec2=generate_sun_array(year,48,meta);
%sunvec2(:,:,2)=mod(sunvec2(:,:,2)+180,360);
zenith2=squeeze(sunvec2(:,11:41,3)');
azim2=squeeze(mod(sunvec2(:,11:41,2)'+180,360));
showplot=0;

if all
    q=1:length(data);
else
    q=p;
end
for k=q
    if preload    
        mornedg=circshift(squeeze(edges(k,:,:,1)),-1,2);
        afteredg=circshift(squeeze(edges(k,:,:,2)),1,2);
    else
        mornedg=squeeze(edges(k,:,:,1));
        afteredg=squeeze(edges(k,:,:,2));
    end
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
    zenith=sunarray(:,150:800,3)';
    azim=sunarray(:,150:800,2)';
    edge=diff(mesh');
    [col,row]=find(edge);
    col=col/1000*48-9;
    simed(k,1:length(row),1)=row';
    simed(k,1:length(col),2)=col';
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



    %x=60:300;
    %Ymorn=polyval(pmorn,x);
    %Yafter=polyval(pafter,x);
    %X=[60, leftedge, leftedge, rightedge, rightedge,300];
    %Y=[leftheight,leftheight,centerheight,centerheight,rightheight,rightheight];
    if showplot
        figure(100);
        imagesc(squeeze(data(k,:,10:40))'); hold on;
        plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
        figure(101);
        imagesc(squeeze(corr(k,:,:))'); hold on;
        plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
        figure(102);

        imagesc((mornedg|afteredg)');hold on;
        plot(row(col<=14),col(col<=14),row(col>14),col(col>14));
        axis([1 days 1 31])
        mornX=zenith2((mornedg|afteredg)');
        mornY=azim2((mornedg|afteredg)');
        figure(104);
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
        figure
        imagesc(mesh');
        %current problem is when only one sided 
        % problem with ouliers
        set(gca, 'YDir','reverse')
    end
end
end
