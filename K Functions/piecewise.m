function out = piecewise(varargin)

    out = nan(size(varargin{1}));
    idx = false(size(varargin{1}));
    conditions = floor(nargin/2);

    for i = 1:2:conditions
        cond = varargin{i};
        val = varargin{i+1};
        out(cond) = val(cond);
        idx = idx | cond;
    end
    other = varargin{end};
    if length(other) == 1
        out(~idx) = other;
    else
        out(~idx) = other(~idx);
    end

end