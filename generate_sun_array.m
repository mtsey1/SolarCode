function time_mat = generate_sun_array(start_year,resolution)
    

%%format inputs for sun_position script
location.latitude = -37.8136;
location.longitude = 144.9631;
location.altitude = 0;
time.UTC = 10;
% time.UTC = 0;
% location.latitude = 47.50685;
% location.longitude = 7.519069;
% location.altitude = 400;



init = datenum(start_year,0,0,0,0,0);

time_mat = zeros(365,resolution,3);

for i=1:365
    for j=0:resolution-1
        time.year = start_year;
        time.month = 1;
        time.day = i;                   % "Jan 32" = "Feb 1"
        time.hour = floor (j * 24 / resolution);
        time.min = 60 * (j * 24 / resolution - time.hour);
        time.sec = 0;

        sun = sun_position(time,location);
        time_mat(i,j+1,2)=sun.azimuth;
        time_mat(i,j+1,3)=sun.zenith;

    end
end




end
