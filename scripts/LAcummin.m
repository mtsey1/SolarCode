% cummin(x,dim)
%
% Explicit formulation of cumulative min, since Matlab 2013b doesn't have it.
%
% Only works for  dim==1,  dim==2  or omitted, and 1- or 2-D  x.

function y=LAcummin(x,dim)
if nargin > 1 || length(size(x)) > 2
    if nargin == 1
        dim = find(size(x) > 1,1);
	if isempty(dim)
	    dim = 1;
	end
    end
    if dim == 1
        y = cummin_(x);
    elseif dim == 2 && length(size(x))==2
        y = cummin_(x')';
    else
	perm = 1:length(size(x));
	perm(dim) = 1;
	perm(1) = dim;
	y = permute(cummin_(permute(x, perm)),perm);
    end
else
    if size(x,1) == 1
	y=cummin_(x')';
    else
	y=cummin_(x);
    end
end
