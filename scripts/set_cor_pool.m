function cor_pool = set_cor_pool (rectangles, valid_days)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Calculate the power consumed by the pump at each time
 % Consider turn on/off within a slot (fractional on/off times) and
 % alternative day-of-week settings.
  if nargin < 2
    valid_days = 1:365;
  end
  cor_pool = zeros(48, valid_days(end))';
  for r = 1:length(rectangles)
    on_off = rectangles(r).on_off;
    st = on_off(1); turn_on = on_off(2); en = on_off(3); turn_off = on_off(4);
    ht = rectangles(r).power;
    missed = rectangles(r).missed;

    days = st:en-1;
    alt_days = [];
    if isfield (rectangles(r), 'alt_days') && sum (rectangles(r).alt_days)
      idx = ismember (mod1 (days, 7), rectangles(r).alt_days);
      alt_days = days(idx);
      days = days(~idx);
    end
    if turn_on < turn_off
      burst = ceil(turn_on):floor(turn_off);        % turn_off must be last
    else                                        % as (1:end-1) is used
      burst = [ceil(turn_on):size(cor_pool,2), 1:floor(turn_off)];    % below
    end
    b = burst(1:end-1);
    d = valid_days (setdiff (days, missed));
    if isempty (alt_days)
      cor_pool(d,b) = min (cor_pool(d,b), -ht);
      frac = ht * (turn_on - floor (turn_on) - 1);
      if frac ~= 0
        if turn_on < 1
          turn_on = turn_on + size(cor_pool, 2);
        end
        cor_pool(d, floor (turn_on)) = min (cor_pool(d, floor (turn_on)), frac);
      end
      frac = ht * (floor (turn_off) - turn_off);
      if frac ~= 0
        if turn_off < 1
          turn_off = turn_off + size(cor_pool, 2);
        end
        cor_pool(d, floor (turn_off)) = min (cor_pool(d, floor (turn_off)), frac);
      end
    else
      % TODO: fix for the case that alt_days is on fewer hours.
      cor_pool(d,b) = min (cor_pool(d,b), -ht);

      turn_on  = rectangles(r).alt_hrs(1);
      turn_off = rectangles(r).alt_hrs(2);
      if turn_on ~= turn_off
        if turn_on < turn_off
          burst = ceil(turn_on):floor(turn_off);
        elseif turn_off > turn_on
          burst = [ceil(turn_on):size(cor_pool,2), 1:floor(turn_off)];
        end
        b = burst(1:end-1);
        b(b==0) = size(cor_pool,2);
        d = valid_days (setdiff (alt_days, missed));
        cor_pool(d, b) = min (cor_pool(d,b), -ht);
      end
    end
  end
end

