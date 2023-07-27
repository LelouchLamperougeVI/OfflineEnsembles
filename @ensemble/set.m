function set(obj,varargin)
if isempty(obj.ops) || nargin<2
    ops.sig = .2;
    ops.e_size = 5;
    ops.thres = 0.75; %threshold for cutting clusters
    
    obj.ops=ops;
    if nargin<2
        return;
    end
end

count = 1;
while(count < length(varargin))
    switch lower(varargin{count})
        case 'sig'
            obj.ops.sig = varargin{count + 1};
        case 'thres'
            obj.ops.thres = varargin{count + 1};
        case 'e_size'
            obj.ops.e_size = varargin{count + 1};
        otherwise
            error(['''' varargin{count} ''' is not a valid parameter']);
    end
    count = count + 2;
end