function [capacity, az, ze, dy, mismatch, funcs] ...
    = findSolarOrientation (seen, sunPos, SamPerDay, ...
                            solar_range, data, s, hot_days, cold_days)

  funcs.angle_coefficient = @angle_coefficient;
  funcs.solar_mismatch = @solar_mismatch;

  seen = seen(:, :, solar_range);
  [mx, mx_pos] = max(squeeze (seen(1,:,:)), [], 2);
  mx(hot_days)  = max (mx(hot_days));
  mx(cold_days) = max (mx(cold_days));

  % Only consider data points at least 25% of peak
  % to exclude both times of high demand and ambient light near dawn/dusk
  big = bsxfun(@gt, squeeze (seen(1,:,:)), mx/4);

  % Find times where the generation is big enough to be trustworthy.
  % This avoids fitting morning/evening (subject to shadows and
  % reflections), and also times with regular daily use, such as the
  % "morning routine" or a pool pump.
  % This is used for finding a lower-bound on capacity,
  % and so it more conservative than the 25% used for the matching.
  big_all = bsxfun (@gt, -data, -0.75 * min(data, [], 2));
  sun_pos.pp = s.pp(big_all);
  sun_pos.s1 = s.s1(big_all);
  sun_pos.s2 = s.s2(big_all);

  % solar_mismatch() will ensure that  cap  is big enough to account
  % for all of the observations recorded in few_data.
  % We could use all timesteps, but for efficiency,
  % we only consider the most negative readings over
  % contiguous 5-day intervals.
  if length(data) > 1
      % Initial estimates, chosen heuristically:
      % Capacity is the maximum observed nett generation
      % Azimuth is roughly proportional to the square of
      % the time of the peak nett generation (with midday being 0),
      % clipped to be between due east and due west
      cap = max(mx);
      az = (((mean(mx_pos) + solar_range(1)-1) * -360/SamPerDay) + 180);
      az = az * abs(az/10);
      az = max (-90, min (90, az));
      ze = 10;		% ~ sun's summer zenith angle in Melbourne
  else                % initialise from previous run
      cap = data;
      az = smr;
      ze = rr;
  end

  options = optimoptions('fmincon', 'Display', 'off', 'Algorithm','active-set');

  max_seen = double (max (seen(:)));
  cost = solar_mismatch (double ([az, ze]), sunPos, ...
                               double (seen), big, max_seen, data(big_all), sun_pos, s.location.latitude);
  if ~isfinite (cost)
    [az, ze, big] = find_feasible (az, ze, sunPos, ...
                                         double (seen), big, ...
                                         max_seen, data(big_all), sun_pos, s.location.latitude);
  end
  old_cost = 1e90;
  
%   for i = -179:180
%       for j = 1:90
%           CostMap(j,i+180)=solar_mismatch([i,j],sunPos, double (seen),big, max_seen, data(big_all), sun_pos, s.location.latitude);
%       end
%   end
%   close all
%   surf(CostMap);
%   save('costmap.mat','CostMap');

  
  
  for n = 1:2
    tic
    % Lower limit of 1 for ze, as ze=0 gives plateau for az.
        [X, cost] = fmincon (@(X) solar_mismatch (X, sunPos, double (seen),...
                                             big, max_seen, data(big_all), sun_pos, s.location.latitude), ...
                         double ([az, ze]), ...
                         [], [], [], [], [-180, 1], [180, 45], [], ...
                         options);
    if isnan(X(1))||isnan(X(2))
        fprintf('X is NaN\n');
        while 1==1
            pause(0.1);
        end
    end
    az  = X(1);
    ze  = X(2);
    [~, ~, cap] = solar_mismatch ([az, ze], sunPos, double (seen), big, ...
                          max_seen, data(big_all), sun_pos, s.location.latitude);

%    fprintf ('pass 1: %g\n', toc);
    if ~isfinite (cap)
      tic
      [X, cost] = fmincon (@(X) solar_mismatch (X, sunPos, double (seen),...
                                               big, max_seen, data(big_all), sun_pos, s.location.latitude), ...
                           double ([az, ze]), ...
                           [], [], [], [], [-90, 1], [90, 45], [], ...
                           options);
      az  = X(1);
      ze  = X(2);
      [~, ~, cap] = solar_mismatch ([az, ze], sunPos, double (seen), big, ...
                            max_seen, data(big_all), sun_pos, s.location.latitude);
      %fprintf ('\t\t\tpass 2: %g\n', toc);
    end
    if cost > 0.99 * old_cost && cost <= old_cost
      break;
    end
    old_cost = cost;
  end
  gen   = cap*angle_coefficient(sunPos, az, ze);

  capacity = cap;
  mx_start = 5;
  dy = squeeze (seen(1, :,:)) - gen;
  mismatch = max(max(abs(dy(:,mx_start:end-mx_start))));
  mismatch = mismatch / max(gen(:));
end


function [cost, dy, cap, gen] = solar_mismatch (X, sunPos, seen, big, ...
                                                max_seen, data, sun_pos, lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%my code%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


  OldTest = 1;
  eps = 0.40905*180/pi;
  if OldTest
     ze  = X(2);
  if 1==1
    az  = X(1);
  else
    t   = X(1);
    if t == 0
        t = 0.001;
    end
    az = t-asind((cosd(ze)/sind(ze))*((sind(t).*cosd(lat-eps))/(sqrt(1-(cosd(t).*cosd(lat-eps)).^2))));
  end
    
  if ~isreal(az)
    cost = inf;
    cap = 0;
    dy = [];
    return
  end 
  
  gen_cap  = angle_coefficient(sunPos, az, ze);
  cap_max = squeeze (seen(1, big(:))) ./ gen_cap(big(:))';
%     fprintf('t is %f\n', t)
%     fprintf('ze is %f\n', ze)
   %s1 is zenith, s2 is azimuth
  capFactor = max (0, cosd (ze) * sun_pos.s1 ...
                     + sind (ze) * cosd (sun_pos.pp - az) .* sun_pos.s2);
  %find the potential output of the panel for all t
  
  
%   t=sun_pos.pp;
%   capFactor = max(0.000001, cosd(t) * cosd(lat-eps) * cosd(ze)  ... 
%       + sqrt(1-(cosd(t) .* cosd(lat-eps)).^2) * sind(ze) .* cosd(t-az));
  
    %fprintf('capFactor is %f\n', capFactor)
    %fprintf('az is %f\n', az)
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%/my code%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%My Method%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if ~OldTest
    az = X(1);
    ze = X(2);
    %fprintf('x = %f, y = %f\n', x,y);
    gen_cap  = angle_coefficient(sunPos, az, ze);
    cap_max = squeeze (seen(1, big(:))) ./ gen_cap(big(:))';
    Xrot = [1,0,0;0,cosd(az),sind(az);0,-sind(az),cosd(az)];
    Yrot = [cosd(ze),0,-sind(ze);0,1,0;sind(ze),0,cosd(ze)];
    %fprintf('Xrot = %f, Yrot = %f\n', Xrot,Yrot);
    Panel_Norm = Xrot*(Yrot*[0;0;1]);

    sunXYZ = [sind(sun_pos.s1).*cosd(sun_pos.pp),cosd(sun_pos.s1).*cosd(sun_pos.pp), sind(sun_pos.pp)];

    capFactor = max(0,((sunXYZ*Panel_Norm)));
    %fprintf(' avg capFactor = %f\n', mean(capFactor));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%/My Method%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  % Find minimum capacity that causes no half-hour to have negative
  % consumption.
  % The factor of 2 comes because cap is in kW,
  % and elements of  data  are in kWh per half hour.
  cap_min = -2 * data ./ capFactor;

  cap_max = 2 * max (cap_max);
<<<<<<< HEAD
  cap_min = max (cap_max / 10, max (cap_min));
=======
  mx = max (cap_max(isfinite (cap_max)));
  if mx < 4
    cap_max = 2 * mx;
  else
    old_mx = mx + 1;
    while old_mx ~= mx
      m = mean (cap_max(cap_max < mx));
      sd = sqrt (var (cap_max(cap_max < mx)));
      old_mx = mx;
      mx = m + 3 * sd;
    end
    % Ignore top two values
    cap_max = 2 * max (cap_max(cap_max < mx));
  end
  cap_min = cap_min(isfinite (cap_min));
  cap_min = max (cap_min(cap_min < mean (cap_min) + 3 * sqrt (var (cap_min))));
  cap_min = min (cap_max / 10, max (cap_min)/2);
>>>>>>> 4919d7713ba20f79c8c0e1e9a7655ebd32b98a8c
  cap_min = double (cap_min);

                                     % sum([]) = 0
  if ~isfinite (az + ze) || ~isfinite (sum (cap_max + cap_min))
    cost = inf;
    cap = 0;
    dy = [];
    return
  end

  if cap_max > cap_min
    cap = cap_min:max (0.01, (cap_max - cap_min) / 100):cap_max + 0.01;
  elseif cap_min > cap_max
    cap = cap_min;          % smallest feasible power
  else
    if ~isempty (cap_max)
      cap = cap_max;    % X:0:Y yields an empty range, so avoid it.
      cap_min = cap;
    else
      cap = cap_min;
      cap_max = cap;
    end
  end

  gen = bsxfun (@times, cap', gen_cap(big)') / 2; % 48 samples
  dy = zeros (size (gen, 1), size (seen(1,:), 2));
  dy(:, big) = bsxfun (@minus, squeeze (seen(1, big)), gen);
  f = (dy > 0);
  if any (f(:))
    value = 100 * (seen ./ max (seen(:))) .^ 2;
    dy(f) = 2 * dy(f) .* value(ceil (find (f) / size (dy, 1)));

    ddy = zeros (size (dy));
    ddy(:, big) = bsxfun (@minus, squeeze (seen(2, big)), gen);
    f = (ddy > 0);
    if any (f(:))
%       dy(f) = dy(f) + 8 * ddy(f);
%       ddy(:, big) = bsxfun (@minus, squeeze (seen(3, bigm)), gen);
%       f = (ddy > 0);
      dy(f) = dy(f) + 900 * ddy(f);
    end
  end

  cost = diag (dy(:,:) * dy(:,:)');
  [~, idx] = min (cost);
  cap  = cap(idx);
  cost = cost(idx);

  cost = cost * max (1.1, cap / max_seen);

  tmp = gen(idx, :);
  gen = zeros (size (big));
  gen(big) = tmp;
end

%%%%%%%%%% non-linear constraint - used for tparam %%%%%%%
function [c,ceq] = nonlcon(x)
    %x(1) is t, x(2) is ze

    if x(2)<abs(27.4165*sind(x(1)))
        c = 1;
        ceq = 1;
        return 
    else
        c = -1;
        ceq = 0;
        return 
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s


%>>>>>>> parent of 41d5bd7... t-parameterisation WORKS - warning: this will mess with other parameterisations
function ac = angle_coefficient (sun, az, ze)
% Fraction of irradience actually received by solar panel,
% given angle of sun and orientation of solar panel.
% Optimzed for  sun  being a vector, and  az and ze being scalars,
% but all/any can be vectors, provided dimensions either match or are 1.
    if numel (sun) == numel ([sun.zenith])
        zenith  = reshape([sun.zenith],  size(sun));	% list -> vector
        azimuth = reshape([sun.azimuth], size(sun));
    else
        zenith  = [sun.zenith];	% list -> vector
        azimuth = [sun.azimuth];
        if size(sun,2) > 1
            zenith  = zenith';
            azimuth = azimuth';
        end
    end

    s1 = cosd(zenith);
    s2 = sind(zenith).*cosd(azimuth - az);

    p1 = cosd(ze);
    p2 = sind(ze);

    ac = max(0, bsxfun(@times, s1, p1) + bsxfun(@times, s2, p2));
end

function [az, ze, big] = find_feasible (az, ze, sunPos, ...
                                              seen, big, ...
                                              max_seen, data, sun_pos, lat)
  cost = solar_mismatch (double ([az, 3]), sunPos, double (seen), ...
                               big, max_seen, data, sun_pos, lat);
  if isfinite (cost)
    ze = 3;
    return;
  end
  cost = solar_mismatch (double ([(az+90)/2, ze]), sunPos, double (seen), ...
                               big, max_seen, data, sun_pos, lat);
  if isfinite (cost)
    az = (az+90)/2;
    return;
  end
  cost = solar_mismatch (double ([(az-90)/2, ze]), sunPos, double (seen), ...
                               big, max_seen, data, sun_pos, lat);
  if isfinite (cost)
    az = (az-90)/2;
    return;
  end

  gen_cap  = angle_coefficient(sunPos, az, ze);
  big = big & gen_cap;
end
