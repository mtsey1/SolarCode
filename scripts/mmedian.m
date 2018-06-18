function m = mmedian (A, dim)
    if isempty(A)
        m = NaN;
    else
	if exist('dim', 'var')
	    m = median(A, dim);
	else
	    m = median(A);
	end
    end
