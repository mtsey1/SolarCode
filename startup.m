ss = get (0, 'ScreenSize');
if isequal (ss, [1 1 1366 768])
  set(0, 'DefaultFigurePosition', [967, 404, 400, 284]);
elseif isequal (ss, [1 1 1536 864])
  set(0, 'DefaultFigurePosition', [1137, 500, 400, 284]);
end
format compact