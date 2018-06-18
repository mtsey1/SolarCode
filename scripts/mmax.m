function m = mmax (A, dim)
    if isempty(A)
        m = -Inf;
    else
	if exist('dim', 'var')
	    m = max(A, [], dim);
	else
	    m = max(A);
	end
    end
