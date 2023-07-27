function [idx, thres] = running(vel)
% Find running epochs

thres = diff(sort(abs(vel)));
thres = round(thres * 2^16) ./ 2^16;
thres = unique(thres);
thres = mean(thres(1:2));

idx = abs(vel) > thres;