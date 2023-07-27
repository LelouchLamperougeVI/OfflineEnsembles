function stacks = make_stacks(behaviour, deconv, varargin)
% Build unit response curves over position and trials
%
% Parameters
%   'bins', 50 (default)
%       number of spatial bins
%
%   'sd', 4 (default)
%       smoothing kernel s.d. in cm


% parse options
ops = parse_input(varargin);

% build behaviour vectors
pos = behaviour.unit_pos;
vel = behaviour.unit_vel;
trials = cat(1, 1, behaviour.trials, length(pos) + 1);
trials = repelem(1:(length(trials) - 1), diff(trials));

% remove immobile epochs and incomplete laps
idx = running(vel);
idx(1:(behaviour.trials(1) - 1)) = false;
idx(behaviour.trials(end):end) = false;
pos = pos(idx);
vel = vel(idx);
trials = trials(idx) - 1;
deconv = deconv(idx, :);

stacks.behaviour.pos = pos;
stacks.behaviour.vel = vel;
stacks.behaviour.deconv = deconv;

% build stacks
edges = linspace(min(pos), max(pos), ops.bins + 1);
pos = discretize(pos, edges);

[idx, ~] = find(deconv > 0);
temp = deconv(deconv > 0);
neur_id = repelem(1:size(deconv, 2), sum(deconv>0, 1));
pos = pos(idx);
trials = trials(idx);

Pt = accumarray([trials(:) pos(:)], 1);
raster = accumarray([trials(:) pos(:) neur_id(:)], temp, [max(trials) max(pos) max(neur_id)]);
stack = squeeze(sum(raster, 1, 'omitnan') ./ sum(Pt, 1));
raster = raster ./ Pt;
raster(isnan(raster) | isinf(raster)) = 0;

% run smoothing
sd = ops.sd / behaviour.vr_length * ops.bins;

temp = reshape(permute(raster, [2 1 3]), [ops.bins, (length(behaviour.trials) - 1) * size(deconv, 2)]);
temp = fast_smooth(temp, sd);
temp = permute(reshape(temp, [ops.bins, length(behaviour.trials)-1, size(deconv,2)]), [2 1 3]);
raster = mat2cell(temp, length(behaviour.trials)-1, ops.bins, ones(1, size(deconv,2)));

stack=fast_smooth(stack, sd);
stack=(stack - min(stack)) ./ range(stack);

% return
stacks.raster = squeeze(raster);
stacks.stack = stack;



function ops = parse_input(inputs)
ops.sd = 4;
ops.bins = 50;

idx = 1;
while(idx < length(inputs))
    switch lower(inputs{idx})
        case 'bins'
            ops.bins = inputs{idx+1};
        case 'sd'
            ops.sd = inputs{idx+1};
        otherwise
            error(['''' inputs{idx} ''' is not a valid parameter']);
    end
    idx = idx + 2;
end