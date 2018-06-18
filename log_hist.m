function [mm, xx] = log_hist (varargin)
  [NN, XX] = hist (varargin{:});

  % convert bars to log scale
  MM(NN>0) = log (NN(NN>0));

  % shift bars so non-zero counts appear positive, and zero counts appear zero
  MM(NN>0) = MM(NN>0) - min(MM) + 1;
  MM(NN==0) = 0;

  bar (XX, MM);
  
  if nargout > 0
    mm = MM;
    xx = XX;
  end
end
