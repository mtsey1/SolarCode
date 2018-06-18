% skipNaN (func, A [, dimension])
% Emulates  func,  but skipping all NaN values in A,
% where  func  is @mean, @sum, @max, @min or similar,
% which takes a single matrix and an optional "dim" argument.
function m = skipNaN (func, A, varargin)
    s = size(A);
    % If no dimension specified, and only one non-singleton, use that.
    if length(varargin) == 0 && length(A(:)) == max(s)
	m = func (A(find(~isnan(A))));
    else
	if length (varargin) > 0
	    dim = varargin{1};
	else
	    dim = find (s > 1, 1);	% Case of all singletons handled above
	end
        if length(s) == 2
	    if dim == 1
		m = zeros(1,s(2));
		for i = 1:s(2)
		    m(i) = func(A(find(~isnan(A(:,i))),i));
		end
	    else
		m = zeros(s(1),1);
		for i = 1:s(1)
		    f = find(~isnan(A(i,:)));
		    if ~isempty(f)
			m(i) = func(A(i,f));
		    else
		        m(i) = NaN;
		    end
		end
	    end
	else
	    % General case, involves reshaping  A, which I assume is expensive.
	    % This may call  func  with an empty matrix -- a problem for median.
	    % To overcome this, call with a shim function that tests for empty.
	    m = vecfun(@(x)(func(x(find(~isnan(x))))), A, dim);
	end
    end
