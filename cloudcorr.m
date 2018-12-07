function correlation=cloudcorr(s,data,meta)
n=60; %length of correlation
width=3; % width of correlation
users=length(s.solar_users);
shadowinglist=[20,146,179,208,233,238,239,240,243,259,266,277,284,292];
correlation=zeros(length(shadowinglist),365,31);
for i=7
    for k=0:364
        for j=0:(31-width)
            [sun,panel]=sunriseset(meta.location.latitude,meta.location.longitude,11,s.solar_az(i),(s.solar_ze(i)),0,s.dark_end);

            % w=data(s.solar_users(i),mod(((0:n)+k),365)+1,(1:3+j+s.dark_end));
            q=(0:(width-1))+j+s.dark_end;
            m=data(s.solar_users(i),mod(((0:n)+k),365)+1,q);
            a=double(squeeze(m));
            % v=s.solar_by_pc(s.postcode(i)==meta.pclist,mod(((0:n)+k),365)+1,1:3+j);
            p=s.solar_by_pc(s.postcode(i)==meta.pclist,mod(((0:n)+k),365)+1,(1:width)+j);
            b=squeeze(p);
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
    cap=s.solar_cap(i);
    post=s.postcode(i);
    correlate=squeeze(correlation(i,:,:));
    corrgrad=gradient(correlate);
    corrgrgrad=gradient(corrgrad);
    datai=squeeze(data(i,:,10:40));
    datagradi=gradient(datai);
    morningshad=(((datai>-cap/10)&(correlate>-0.3)&(corrgrad<0)&(corrgrgrad<0))');
    afternoonshad=(((datai>-cap/10)&(correlate>-0.3)&(corrgrad>0)&(corrgrgrad<0))'); 
    figure(100)
    image(squeeze(correlation(i,:,:))','CDataMapping','scaled')
    hold on
    sunriseset(meta.location.latitude,meta.location.longitude,11,s.solar_az(i),(s.solar_ze(i)),1,s.dark_end);
    axis([1 365 1 31])
    figure(101); 
    imagesc((datai>0)')
    %imagesc(imgaussfilt(squeeze(data(i,:,10:41))',[2 0.5]));
    figure(102);
    imagesc(afternoonshad)
    hold on
    sunriseset(meta.location.latitude,meta.location.longitude,11,s.solar_az(i),(s.solar_ze(i)),1,s.dark_end);
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
    sunriseset(meta.location.latitude,meta.location.longitude,11,s.solar_az(i),(s.solar_ze(i)),1,s.dark_end);
    axis([1 365 1 31])
    %count(i)=nnz(squeeze((correlation(i,:,:)>-0.2)&(correlation(i,:,:)~=0))')
    close([100 104 102]);
end
        
end

% maybe look at corr2 for some range of hours instead of the vectors used 
