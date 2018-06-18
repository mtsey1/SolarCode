function best = find_rect (d)
% BEST = FIND_RECT (D)
% Return a triple (x_min, x_max, area) such that a rectangle of size
% area spanning elements x_min to x_max of d, entirely above the graph
% of d and below the lower of the two peaks.
    p = min(d([1,end]));
    d = d(2:end-1);
    dd = sort(d);
    best = [0, 0, 0];
    for j = length(dd):-1:1
        r = ranges(find(d <= dd(j)));
        [len, pos] = max(r(2,:)-r(1,:)+1);
        height = p - dd(j);
        area = len * height;
        if area > best(3)
            best = [r(:,pos)', area];
        else
            if len * p < best(3)
                break
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function rect = find_rect (peak, d, idx, inv)
%    rect = find_me;
%
%    rect1 = find_rect (peak, d(something), idx(something), inv(something));
%    if rect1(3) > rect(3)
%        rect = rect1;
%    end
%
%    rect1 = find_rect (peak, d(something), idx(something), inv(something));
%    if rect1(3) > rect(3)
%        rect = rect1;
%    end
%
%end


%    primPk = zeros(1,size(data,1));
%    for k = 1:size(data,1)
%        if data(k,2) < data(k,1)
%	    primPk(k) = 1;
%	else
%	    primPk(k) = find(diff(data(k,:)) < 0, 1);
%	end
%    end
%    cm = cummin(data,2);
%    [tmp p] = max(data - cm, [], 2);
%    for k = find(p <= primPk)
%        [tmp pp] = max(data(k,primPk(k):end) - cummin(data(k,primPk(k):end)), [], 2);
%	p(k) = primPk-1 + pp;
%    end
%end
