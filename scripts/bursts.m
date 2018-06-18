function b = bursts (cv)
% Record properties of bursts of above-average power, of various lengths.
% Record Power, duration, time-of-year
  big_min = cv;
  all_bursts = [];
  for i = 1:max_burst
    if i ~= 1
      big_min = min(big_min(1:end-1), cv(2:end));
    end
    pre = big_min(1:end-2);
    post = big_min(3:end);
    height = big_min(2:end-1) - max (pre, post);

    times = find (height > abs (pre - post));
    new_bursts = [ones(length(new_burst),1), height(times), mod(times+1, size(cv,1)), ];
    all_burst = [all_bursts; new_bursts];
  end


  % Cluster the bursts.  How?
end
