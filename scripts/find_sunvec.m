function [errorvec, errorsize]=find_sunvec(vec,showplot)
year = 2013;
pos=find(vec);
for q=1:6
pos(end+1)=NaN;
end
thresh1=1.5;
thresh2=0.5;
sunvec=generate_sun_array(year,48);
[yy,xx]=ind2sub(size(vec),pos);
errorvec=zeros(size(vec));
i=1;
fprintf('%d \n',length(pos))
while(~isnan(pos(i)))
    %another way may be to find the loccation with the minimum distance
    %away from it 
    %           Ze                                                  Az
    dist=(abs(sunvec(:,11:41,2)-sunvec(xx(i),yy(i)+10,2)))+5.*(abs(sunvec(:,11:41,3)-sunvec(xx(i),yy(i)+10,3)));
    [res, index]=mink(reshape(dist,[],1),7); 
    logicsunvec=zeros(365,31);
    pos2=index;
    for k=1:5
        logicsunvec(ind2sub(size(dist),index(k)))=1;
    end
    %{
    logicsunvec=(abs(sunvec(:,11:41,2)-sunvec(xx(i),yy(i)+10,2))<thresh1).*(abs(sunvec(:,11:41,3)-sunvec(xx(i),yy(i)+10,3))<thresh2);
    pos2=find(logicsunvec);
    %}
    sum=0;
    for j=1:length(pos2)
        [y,x]=ind2sub(size(vec'),pos2(j));
        if vec(x,y)
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
    end   
    %}
    i=i+1;
    if showplot
        figure(50);
        imagesc(vec)
        figure(51);
        imagesc(logicsunvec')
        figure(52);
        imagesc(errorvec)
    end
    
end
fprintf('%d \n',length(pos))
errorsize=length(find(errorvec))/length(pos);
figure;
imagesc((vec|errorvec)+errorvec);
end
