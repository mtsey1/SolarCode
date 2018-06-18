% Copyright (C) 2016 Lachlan Andrew

% -*- texinfo -*- 
% @deftypefn {} {@var{retval} =} fridge (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: Lachlan Andrew <lachlan.andrew@monash.edu>
% Created: 2016-12-13

function [model, err] = fridge (measured, amp, phase, freq, duty)
  if nargin == 1
    measured = double (measured);
    wrap = [measured(2); measured(:); measured(end-1)];
    means = (wrap(1:end-2) + wrap(3:end)) / 2;
    %diffs = abs (wrap(1:end-2,:) - wrap(3:end, :));
    %peaks = wrap(2:end-1,:) > means + 2 * diffs;
    peaks = [false; wrap(2:end-1)>means; false];
    
    maxes = max (wrap(1:end-2), wrap(3:end));
    mins  = min (wrap(1:end-2), wrap(3:end));
    % TODO: should be <= on one side and < on the other
    %       peaks should have weak inequality on opposite side of troughs
    peakpeaks = [false; wrap(2:end-1)>maxes; false];
    troughs   = [false; wrap(2:end-1)<mins;  false];
    
%{
cv = squeeze (data(i,:,:))';
aa = cv(1:10, 1:4);
plot (aa(:))
aa = aa(:);
ma = mean (aa);
above = find (aa > ma);
below = find (aa < ma);
plot (above, aa(above), below, aa(below), 1:length(aa), aa)
%}
    
    % TODO: filter: Use max of neighbours if centre is much larger
    
    % TODO: Identify "base" samples as those either side of a spike
    % that are almost equal, say to within 10% of the height of the spike.
    
    % mean of samples separated by 3
    mean_3 = (wrap(4:end) + wrap(1:end-3)) / 2;
    % smear left and smear right
    sl = [false; (peaks(3:end-1) & wrap(2:end-2) > wrap(1:end-3) ...
                      &   abs (wrap(2:end-2) - wrap(4:end)) ...
                        > abs (wrap(1:end-3) - wrap(4:end)))];
    sr = peaks(2:end-2) & wrap(3:end-1) > wrap(4:end) ...
                        &   abs (wrap(3:end-1) - wrap(1:end-3)) ...
                          > abs (wrap(4:end)   - wrap(1:end-3));

    % TODO: For those not marked as smeared,
    % but for which the point next to the lower is a "base" sample,
    % mark as smeared also.
    
    % TODO: From fraction of smeared overnight,
    %       estimate the duration of "ON".
    %       (Work out what to do if duration > 1)
    
    % TODO: estimate on/off times of smeared.
    %       interpolate to find on/off times of non-smeared.
    %       if these lie outside peak interval,  do something...
    
    peaks = double (peaks);
    peaks(logical (peaks)) = find(peaks);   % logical -> spaced indices
    peaks(sl) = peaks(sl) - 1 ...
      +  (wrap([false; sl]) - mean_3(sl)) ...
      ./ (wrap([false; sl]) + wrap(sl) - 2*mean_3(sl));
    peaks(sr) = peaks(sr) ...
      +  (wrap(sr) - mean_3([false; sr])) ...
      ./ (wrap(sr) + wrap([false; sr]) - 2*mean_3([false; sr]));
    peaks = max (0, reshape (peaks(2:end-1), size (measured)) - 1);
    
    %p = peaks(1:11, :);
    %gaps = diff (find (p(:)));
    %mg = mean (gaps(gaps >= 2 & gaps <= 4));

    for i = 1:size (measured, 2)
      ons = mod (peaks (peaks(1:11, i) ~= 0, i), size (peaks, 1));
      meas = measured(1:11,i);
      offs = ons + 1;
      amp = 0.05;
      len = length (ons);
      if len == 0
        continue;
      elseif len > 1
        offs(1:end-1) = min(offs(1:end-1), (ons(1:end-1) + ons(2:end)) / 2);
      end
      obj = @(x)(mismatch (meas, x(end), x(1:len), x(len+1:2*len)));
      model = fmincon (obj, [ons; offs; amp], [], [], [], [], [-100*ones(2*len, 1); 0], [length(meas)*ones(2*len, 1); mean(measured(:))]);
      amps(i) = model(end);
      models{i} = model(1:end-1);
      ons = model(1:len);
      offs = model(len+1:2*len);
      amp = model(end);
      [~, gen] = mismatch (meas, amp, ons, offs);
      [~] = abs (gen);  % avoid 'gen unused' warning
      keyboard
    end
    ons = model(1:len);
    offs = model(len+1:2*len);
    amp = model(end);
    %[~, gen] = mismatch (measured, amp, ons, offs);
    return
  end
    
  ons = phase(:);
  offs = freq(:);
  len = length (ons);
  measured = double (measured);
  obj = @(x)(mismatch (measured, x(end), x(1:len), x(len+1:2*len)));
  model = fmincon (obj, [ons; offs; amp], [], [], [], [], [-100*ones(2*len, 1); 0], [length(measured)*ones(2*len, 1); mean(measured)]);
  ons = model(1:len);
  offs = model(len+1:2*len);
  amp = model(end);
  [~, gen] = mismatch (measured, amp, ons, offs);
  return

  raw = zeros (size (measured) .* [300, 1]);
  x = cos (((1:size (raw, 1)) - phase * 300) * 2 * pi * freq / 300);
  x = (x > cos (duty * pi));
  model = sum (reshape (x, [300, size(measured)]));
  figure(3); plot ((1:length (x))/300, x, [1:size(measured, 1); 1:size(measured, 1)], repmat ([1; 0], 1, size(measured,1)), 'r');
  figure(2); plot (model);
  err = sum (model - measured);
end

function [model, err] = fit_am (measured, cycles)
  ons(1) = size (measured, 1) / cycles;
end

function [cost, gen] = mismatch (measured, amp, ons, offs)
  % Generate the data we should have measured, given amplitude and on/off
  % times.
  
  % First, set time slots with no on/off events
  fon  = floor (ons);
  step_on = zeros (size (measured));
  step_on(fon(fon >= 1)) = 1;
  s_on = cumsum (step_on);

  foff = floor (offs);
  step_off = zeros (size (measured));
  step_off(foff(foff >= 1)) = 1;
  s_off = cumsum (step_off);

  gen = s_on - s_off;
  gen = gen - min (gen);    % if {-1, 0}, make it {0, 1}

  % Second, set time slots with a single on/off event
  a = [ons, zeros(size (ons)); offs, ones(size (offs))];
  a = a(a(:,1) > 0 & a(:,1) <= length (measured), :);
  [~, idx] = sort (a(:,1));
  a = a(idx, :);
  d = diff (a(:, 2));
  if any (d == 0)
    cost = 1e30 * sum (d == 0);
    fprintf ('out of bounds\n');
    return;
  end
  
  c = ceil  (a(:,1));
  f = floor (a(:,1));
  d = diff (f);
  if isempty (d)
    isolated = (f >= 1);
  else
    isolated = [d(1)>0; d(1:end-1)&d(2:end); d(end)>0];
    isolated = isolated & (f >= 1);
  end
  iso_on = isolated & a(:,2) == 0;
  gen(f(iso_on)) = 1 + f(iso_on) - a(iso_on, 1);
  iso_off = isolated & a(:,2) == 1;
  gen(f(iso_off)) = a(iso_off, 1) - f(iso_off);
  
  % Finally, set time slots with multiple on/off events
  if length (isolated) > 1
    isolated = [true; isolated];
    bunches = find (isolated(2:end) < isolated(1:end-1));
    if ~isempty (bunches)
      for bunch = bunches'
        i = bunch;
        g = 0;
        while ~isolated(i) && i <= size (a, 1)
          if a(i,2) == 0
            g = g - a(i, 1);    % Earlier on  -> more energy
          else
            g = g + a(i, 1);    % Earlier off -> less energy
          end
          i = i + 1;
          if f(i) > f(i-1)
            gen(f(bunch)) = mod (g, 1);
            g = 0;
            bunch = i;
          end
        end
        %gen(f(bunch)) = mod (g, 1);
      end
    end
  end
  
  gen = gen * amp;
  
  if nargout > 1
    figure(1);
    plot (1:length (measured), measured, 1:length (gen), gen);

    figure(2);
    plot (measured - gen);
    %keyboard
  end
  
  err = measured - gen;
  variability = sum (abs (diff (diff (ons)))) ...
              + sum (abs (diff (offs - ons)));
  cost = sum (abs (diff (err))) + 0.1 * variability;
  fprintf ('%g\n', cost);
end

% fridge (a(:, 1), 0.05, 4.45, 0.455, 0.25)
% fridge (a(:, 2), 0.05, 0.78, 0.41, 0.3);
% fridge (a(:, 3), 0.05, -1.265, 0.481, 0.1);
% fridge (a(:, 4), 0.05, -0.2, 0.4835, 0.2);
% fridge (a(:, 5), 0.05, 0.52, 0.376, 0.2);
% fridge (a(:, 6), 0.05, -0.4, 0.35, 0.3);