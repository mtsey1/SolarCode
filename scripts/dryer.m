function cor = dryer (user_data)
% Identify possible turn-on times for a clothes dryer
% Currently only tested on customer 212, and it doesn't even work for that.
  pwr = 0.6;
  duration = 4;
  dur2 = floor (duration / 2);
  
  smoothed = medfilt1 (user_data, 3, [], 2);
  
  % Find points exceeding surrounds by pwr.
  surround = min ([nan(1, dur2), smoothed(1:end - dur2)], ...
                  [smoothed(floor (duration/2) + 1:end), nan(1, dur2)]);
  surround = reshape (surround, size (user_data));
  candidates = (user_data - surround > pwr ...
              & user_data - surround < 15 * pwr);
  
  % TODO: Recalculate pwr and duration by clustering candidates

  % Isolate candidates
  candidates(1) = (candidates(1) && ~candidates(2));
  candidates(2:end-1) = (candidates(2:end-1) & ~candidates(1:end-2));
  
  %idle = medfilt1 (user_data, 13, [], 2) < 0.1;
  %candidates = candidates & idle;
  
  %candidates (:, 60:end) = false;
  candidates([1:15, 35:end], :) = false;
  
  starts = find (candidates(:));
  overlaps = zeros (duration + 2, length (starts));
  for i = -1:duration+2
    overlaps(i+2, :) = user_data (starts + i - 1);
  end
  figure(107); plot (overlaps);
  
  % Check (usually) not more than once per day
  
  % If less than once per week, correlation with rain  
  
  cor = reshape (candidates, size (user_data));
end