function a = plot_runs (runs, meta)
    tmp = zeros(meta.SamPerDay, meta.Days);
    for j = 1:size (runs, 1)
      if runs(j,3) ~= -1
        tmp(runs(j,3), runs(j,1)+1:runs(j,2)) = runs(j,7) * sign (runs(j,5));
      end
    end
    imagesc (tmp);
end
