function show_heatmap (var, rectangles, cv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  values = zeros (size (cv));
  on_off = [rectangles.on_off];
  for r = 1:length (rectangles)
    values (rectangles(r).burst, on_off(1,r):on_off(3,r)-1) = -var(r);
  end
  imagesc (values);
end

