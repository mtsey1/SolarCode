function [r, str] = ranges(fd)
 % Detect runs of consecutive integers (or reals spaced by 1) in the input
 % The output is a 2-by-k matrix in which the first row is the first element
 % in each of the  k  runs, and the second row is the last element in each run.
    if isempty(fd)
        r = [];
        str = '';
    else
        fd = fd(:)';            % force to row, regardless of shape of fd
        dfd = diff(fd);
        ends = find(dfd ~= 1);  % report [1 2 2 3] as two ranges: 1-2, 2-3
        ends = ends(:)';
        r = [fd([1,ends+1]); fd([ends,end])];
        if nargout >= 2
            str = sprintf('%d-%d ', r);
        end
    end
end
