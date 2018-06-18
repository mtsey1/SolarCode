function show_aircon (aircon, background, colour, weekends)
  if nargin < 3
    colour = [0, 0, 0];
  end
  if nargin < 4
    weekends = [];
  end

    % Display power 'heatmap', preserving any zoom
  xx = get (gca, 'XLim');
  yy = get (gca, 'YLim');
  imagesc (background);
  if all (xx ~= 0) && all (yy ~= 0)
    set (gca, 'XLim', xx);
    set (gca, 'YLim', yy);
  end

  if ~isempty (weekends)
    for i = weekends(:)'
      line ([i, i], [1, size(background, 2)], 'Color', [1, 0.5, 0]);
    end
  end

  if ~isempty (colour)
    for r = ranges (find (aircon))
      [y, x] = ind2sub (size (background), r);
      if x(1) == x(2)
        line(x, y(:) + [-0.5; 0.5], 'Color', colour);
      else
        line ([x(1), x(1)], [y(1)-0.5, size(background, 1)+0.5], 'Color', colour);
        line ([x(2), x(2)], [0.5, y(2)+0.5], 'Color', colour);
        for i = x(1)+1 : x(2)-1
          line ([i, i], [0.5, size(background, 1)+0.5], 'Color', colour);
        end
      end
    end
  end
end
