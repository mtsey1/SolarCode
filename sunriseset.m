function [solar,panel]=sunriseset(lat,long,UTCoff,normalAz,normalZe,plotout,offset)
%takes in the latitude longtitude UT offset and the normalAZ(angle of solar
%panel from facing the equator and normalZe(the angle of the panel above
%the horizontal and whether to output a plot. 
%the output is the sunrise and sunset times for both the sun and the solar
%panel 

solar=zeros(2,365);
longCorr = 4*(long - 15*UTCoff);                         % longitudinal correction
days = 1:365;
B = 360*(days - 81)/365;
EoTCorr = 9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);  % Equation of Time correction
solarCorr = longCorr + EoTCorr;

delta = 23.45*sind(360*(days + 287)/365);    % Solar declination


solar(1,:) = 2*(12 - acosd(-tand(lat)*tand(delta))/15 - solarCorr/60);
solar(2,:)  = 2*(12 + acosd(-tand(lat)*tand(delta))/15 - solarCorr/60);
%%%%---Method 1------%%%%
%may need to transform panel zenith by 90 deg
%my need to try negative of the panel Az
if (lat<0)
    Az= 180-normalAz;
end   
a=cosd(lat)/(sind(Az)*tand(normalZe))+sind(lat)/tand(Az);
b=tand(delta)*(cosd(lat)/tand(Az)-sind(lat)/(sind(Az)*tand(normalZe)));
k1=(a*b-sqrt(a^2-b.^2+1))/(a^2+1);
k2=(a*b+sqrt(a^2-b.^2+1))/(a^2+1);
ws1=acosd(k1);
ws2=acosd(k2);
if (normalAz<0)
    panel(1,:)=2*(12-ws2/15-solarCorr/60);
    panel(2,:)=2*(12+ws1/15-solarCorr/60); 
else
    panel(1,:)=2*(12-ws1/15-solarCorr/60);
    panel(2,:)=2*(12+ws2/15-solarCorr/60); 
end

if plotout
    plot(days, solar(1,:)-offset, days, solar(2,:)-offset, 'LineWidth', 2)
    hold on
    plot(days, panel(1,:)-offset, days, panel(2,:)-offset, 'LineWidth', 2)
    title('Sunrise and Sunset')
    xlabel('Day of Year')
    ylabel('Time')
    legend('Sunrise', 'Sunset','Panel Sunrise','Panel Sunset')
    yticks([0,5,10,15,20,25,30])
    yticklabels({'5:00','7:30','10:00','12:30','15:00','17:30','20:00'})
end

%%%%------Method 2------%%%%
%{
long=long-normalAz;
lat=lat+normalZe;
panel=zeros(2,365);
longCorr = 4*(long - 15*UTCoff);                         % longitudinal correction
days = 1:365;
B = 360*(days - 81)/365;
EoTCorr = 9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);  % Equation of Time correction
solarCorr = longCorr + EoTCorr;

delta = asind(sind(23.45)*sind(360*(days - 81)/365));    % Solar declination

panel(1,:) = 12 - acosd(-tand(lat)*tand(delta))/15 - solarCorr/60;
panel(2,:)  = 12 + acosd(-tand(lat)*tand(delta))/15 - solarCorr/60;

plot(days, panel(1,:), days, panel(2,:), 'LineWidth', 2)
axis([1 365 0 24])
%}
%%%%-----Method 3----%%%%
%{
D=sind(delta)*sind(lat)*cosd(90-normalZe)-sind(delta)*cosd(lat)*sind(90-normalZe)*cosd(normalAz);
A=cosd(delta)*cosd(lat)*cosd(90-normalZe);
B=cosd(delta)*sind(lat)*sind(90-normalZe)*cosd(normalAz);
C=cosd(delta)*sind(90-normalZe)*sind(normalAz);
den=(A.^2+2*A.*B+B.^2+C.^2);
sq=sqrt(A.^2.*C.^2+2*A.*B.*C.^2+B.^2.*C.^2+C.^4-C.^2.*D.^2);
n1=((-A.*D-B.*D-sq)/den);
n2=((-A.*D-B.*D+sq)/den);
d1=(1./C).*(-D+A.^2.*D./den+2*A.*B.*D./den+B.^2.*D./den+A.*sq./den+B.*sq./den);
w1=atan2(n1,d1)*360/(2*pi);
w2=atan2(n2,d1)*360/(2*pi);
panel(1,:)=(12+w1/15-solarCorr/60);
panel(2,:)=(12+w2/15-solarCorr/60);
plot(days, panel(1,:), days, panel(2,:), 'LineWidth', 2)
%}



