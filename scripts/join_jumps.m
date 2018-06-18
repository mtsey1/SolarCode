function [jumps, times, signed_jumps] = join_jumps(me_year)
  % MY_JUMPS -- estimate times and sizes of large jumps
  me_year = me_year(:);
  times = (1:length (me_year) - 1)';
  signed_jumps = me_year(2:end) - me_year(1:end-1);
  jumps = abs (signed_jumps);
  mj = mean (jumps);
  sj = sign (signed_jumps);
  consecutive = (sj(1:end-1) == sj(2:end)) ...
                & ((jumps(1:end-1) + jumps(2:end)) > 2 * mj);
  % Two consecutive slots may be marked "consecutive", so we can't
  % just sum them to get the jump.
    
  % Record start/end of consecutive sequences,
  % *before* marking the last of a sequence as part of the run,
  % so that contiguous runs don't merge.
  first = consecutive & [true; ~(consecutive(1:end-1))];
  last  = consecutive & [~consecutive(2:end); true];
  pairs = find (first & last);
  consecutive (pairs) = false;

  % Flag last of a consecutive run as a "consecutive" too
  consecutive([false; consecutive(1:end-1)]) = true;
  last = [false; last(2:end)];
  
  mid = jumps(consecutive);
  pre = jumps(consecutive(2:end));
  if consecutive(1)
    pre = [0; pre];
  end
  post = jumps([false; consecutive(1:end-1)]);
  if consecutive(end)
    post = [post; 0];
  end
  
  % Handle runs of consecutive signs.
  first = first(consecutive);
  last = [false; last(consecutive)];
  if length (last) > length (first)
    last = last(1:end-1);
  end
  
  take_post    = mid >  post & ((post >= pre) | first) & ~last;
  give_to_post = mid <= post & ((post >= pre) | first) & ~last;
  take_pre     = mid >= pre  & ((pre >= post) | last)  & ~first;
  give_to_pre  = mid <  pre  & ((pre >= post) | last)  & ~first;
  
  rev = find(consecutive);    % map back to index in full matrix
  
  take_post = rev(take_post(1:end-1) & give_to_pre( 2:end));
  take_pre  = rev(take_pre( 2:end)   & give_to_post(1:end-1)) + 1;
  take_post = [pairs; take_post];
  
  times(take_post) = (times(take_post+1) .* jumps(take_post+1) ...
                      + times(take_post) .* jumps(take_post)) ...
                     ./ (jumps(take_post) + jumps(take_post+1));
  jumps(take_post) = jumps(take_post) + jumps(take_post+1);
  signed_jumps(take_post) = signed_jumps(take_post) + signed_jumps(take_post+1);
  times(take_post+1) = 0;
  jumps(take_post+1) = 0;
  %jumps([false, take_post) = 0;
  %signed_jumps([false, take_post]) = 0;
  
  times(take_pre) = (times(take_pre-1) .* jumps(take_pre-1) ...
                     + times(take_pre) .* jumps(take_pre)) ...
                    ./ (jumps(take_pre) + jumps(take_pre-1));
  jumps(take_pre) = jumps(take_pre) + jumps(take_pre-1);
  signed_jumps(take_pre) = signed_jumps(take_pre) + signed_jumps(take_pre-1);
  times(take_pre-1) = 0;
  jumps(take_pre-1) = 0;
  %jumps(take_pre(2:end) = 0;
  %signed_jumps(take_pre(2:end)) = 0;
  
  keep = (times ~= 0) & (jumps ~= 0);
  jumps = jumps(keep);
  signed_jumps = signed_jumps(keep);
  times = times(keep);
  
%   merge = zeros (signed_jumps);
%   merge(consecutive) = signed_jumps(consecutive) + signed_jumps([false, consecutive]);
%   choose_before = merge(1:end-1) > merge(2:end);
%   choose_after = 1;
%   i_want_to_absorb_R = consecutive & (jumps(1:end-1) > jumps(2:end)) ...
%                         & (~consecutive([2:end, end]) ...
%                            | jumps(2:end) > jumps([1,1:end-2]));
end