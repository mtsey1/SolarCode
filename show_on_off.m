%% Copyright (C) 2016 Lachlan Andrew
%% 
%% -*- texinfo -*- 
%% @deftypefn {} {} show_on_off (@var{on_off}, @var{background}, @var{meta})
%%
%% @seealso{}
%% @end deftypefn

%% Author: Lachlan Andrew <lachlan@la-dell>
%% Created: 2016-07-12

function a =show_on_off (oo, cv, meta)
  imagesc (cv);
  
  if isempty (oo)
    return
  end
  oo([1,3],:) = oo([1,3],:) - 0.5;

%  if isfield(rr, 'colour')
%    colour = { rr.colour }';
%    a = cellfun (@isempty, colour);
%    colour(a,:) = repmat ({[1 1 1]'}, 1, sum (a));
%    set (gca, 'colorOrder', [colour{:}]');
%  else
%    set (gca, 'colorOrder', [1 1 1]);
%  end

  s = (oo(4,:) > oo(2,:));   % "simple" case
  if any (s)
    x = [oo(1,s); oo(1,s); oo(3,s); oo(3,s); oo(1,s)];
    y = [oo(2,s); oo(4,s); oo(4,s); oo(2,s); oo(2,s)]-0.5;
    line (x, y, 'Color', 'w');
  end

  s = ~s;   % "split over days" case
  if any (s)
    one = ones (1, sum (s));
    x = [oo(1,s); oo(1,s); oo(3,s); oo(3,s)];
    y1 = [(meta.SamPerDay+1)*one; oo(2,s); oo(2,s); (meta.SamPerDay+1)*one];
    y2 = [one; oo(4,s); oo(4,s); one];
    line ([x, x+1], [y1, y2] - 0.5, 'Color', 'w');
  end
  
  a = 0;
end
