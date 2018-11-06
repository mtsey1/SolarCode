function ave = moving_ave (data, length, dim)
  if numel (data) == 0
    ave = data;
    return
  end
  
  if nargin < 2
    error ('moving_ave requires at least two arguments');
  end
  
  if nargin < 3
    [~, dim] = find (size (data) > 1, 1);
    if isempty (dim)
      dim = 1;
    end
  end
  if length > size (data, dim)
    length = size (data, dim);
  end
  
  cs = cumsum (data, dim);
  ave = zeros (size (data));
  left = floor (length / 2);
  if dim == 1
    ave(left+2:end-left, :) = (cs(length+1:end, :) - cs(1:end-length, :)) / length;
    ave(1:left+1, :) = bsxfun (@rdivide, cs(left+1:length, :), (left+1:length)');
    ave(end-left+1:end, :) = bsxfun (@rdivide, cs(end, :) - cs(end-length+1:end-left-1, :), (length-1:-1:left+1)');
  else
    sz = size (ave);
    ave = reshape (ave, prod (sz(1:dim-1)), sz(dim), []);
    cs  = reshape (cs, size (ave));
    
    ave(:, left+2:end-left, :) = (cs(:, length+1:end, :) - cs(:, 1:end-length, :)) / length;
    ave(:, 1:left+1, :) = bsxfun (@rdivide, cs(:, left+1:length, :), (left+1:length));
    ave(:, end-left+1:end, :) = bsxfun (@rdivide, cs(:, end, :) - cs(:, end-length+1:end-left-1, :), (length-1:-1:left+1));
    
    ave = reshape (ave, size (data));
  end
  
end