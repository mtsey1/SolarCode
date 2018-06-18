% squeezeNaN (func, A [, dimension])
% Emulates  func,  but skipping all NaN values in A,
% where  func  is @mean, @sum, @max, @min or similar,
% which takes a single matrix and an optional "dim" argument.
% Finally "squeezes" out the dimension over which the action was taken,
% but unlike Matlab's "squeeze", does not remove earlier singleton dimensions.
function m = squeezeNaN (func, A, varargin)
    s = size(A);
    % If no dimension specified, and only one non-singleton, use that.
    if length(varargin) == 0 && length(A(:)) == max(s)
	m = func (A(~isnan(A)));
    else
	if ~isempty (varargin)
	    dim = varargin{1};
	else
	    dim = find (s > 1, 1);	% Case of all singletons handled above
	end
        if length(s) == 2
	    if dim == 1
		m = zeros(1,s(2));
		for i = 1:s(2)
		    m(i) = func(A(~isnan(A(:,i)),i));
		end
	    else
		m = zeros(s(1),1);
		for i = 1:s(1)
		    f = (~isnan(A(i,:)));
		    if sum(f)>0
			m(i) = func(A(i,f));
		    else
		        m(i) = NaN;
		    end
		end
	    end
	else
	    % General case, involves reshaping  A, which I assume is expensive.
	    % This may call  func  with an empty matrix -- a problem for median.
	    % To overcome this, the caller can call us with a shim function
	    % that tests for empty.
	    m = vecfun(@(x)(func(x(~isnan(x)))), A, dim);
	end
    end
    if length(s) > 2
	m = reshape(m, [s(1:dim-1), s(dim+1:end)]);
    end
