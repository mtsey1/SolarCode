cc = max (0, cv + (cor(1:size(cv,2),:))');
grid = 100;
a = (0:grid) .^ 2 * max (cc(:)) / (grid .^ 2);
%h = zeros (grid+1, grid);
v = zeros (1, grid);
for z = 1:grid
  first = (cc(1:end-1) > a(z) & cc(1:end-1) < a(z+1));
%  h(:, z) = histc (cc([false, first]), a);
  v(z) = var (cc ([false, first]));
end
%h = bsxfun (@rdivide, h, sum (h));
%figure(1); imagesc (log (h));
v = (v + medfilt1 (v)) / 2;
v = (v + medfilt1 (v)) / 2;
%figure(2); plot (a(2:end), v + medfilt1 (v));
%figure(3); plot (a(2:end), sqrt (v + medfilt1 (v)));