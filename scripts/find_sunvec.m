function [errorvec, errorsize]=find_sunvec(shading,vec,showplot,meta)
% this is crude and currently unused 
year = meta.Year;
pos=find(vec);
temp=length(pos);
showplot=0;
for q=1:12
pos(end+1)=NaN;
end
thresh1=1.5;
thresh2=0.5;
sunvec=generate_sun_array(year,48,meta);
[yy,xx]=ind2sub(size(vec),pos);
errorvec=zeros(size(vec));
i=1;
if showplot
    fprintf('%d \n',length(pos))
end
count=0;
while(~isnan(pos(i)))
    %another way may be to find the loccation with the minimum distance
    %away from it 
    %           Ze                                                  Az
    dist=(abs(sunvec(:,11:41,2)-sunvec(xx(i),yy(i)+10,2)))+5.*(abs(sunvec(:,11:41,3)-sunvec(xx(i),yy(i)+10,3)));
    [res, index]=mink(reshape(dist,[],1),7); 
    logicsunvec=zeros(365,31);
    pos2=index;
    for k=1:7
        logicsunvec(ind2sub(size(dist),index(k)))=1;
    end
    %{
    logicsunvec=(abs(sunvec(:,11:41,2)-sunvec(xx(i),yy(i)+10,2))<thresh1).*(abs(sunvec(:,11:41,3)-sunvec(xx(i),yy(i)+10,3))<thresh2);
    pos2=find(logicsunvec);
    %}
    sum=0;

    for j=1:length(pos2)
        [y,x]=ind2sub(size(vec'),pos2(j));
        if shading(x,y)
            %should only remove once it it has been gone over for many
            %locations
             pos(pos==pos2(j))=[];
        else
            sum=sum+1;
            errorvec(x,y)=1;
        end
    end
    if 2*sum>length(pos2)
        errorvec(yy(i),xx(i))=1;
    else
        count=count+1;
    end
    %}
    i=i+1;
    if 1==0
        figure(50);
        imagesc(vec)
        figure(51);
        imagesc(logicsunvec')
        figure(52);
        imagesc(errorvec)
    end
    
end

errorsize=count/temp;
if 1==0
    fprintf('%d \n',length(pos))
    figure;
    imagesc((vec|errorvec)+errorvec);
end
end
