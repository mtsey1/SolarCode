function v = vvar (A, dim)
    if isempty(A)
        v = NaN;
    else
	if exist('dim', 'var')
	    v = var(A, dim);
	else
	    v = var(A);
	end
    end
