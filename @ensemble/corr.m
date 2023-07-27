function corr(obj)
% build correlation matrix

deconv = obj.twop.deconv;
deconv = fast_smooth(deconv,obj.ops.sig * obj.twop.fs);

deconv(all(isnan(deconv), 2), :) = [];
deconv=(deconv - mean(deconv)) ./ std(deconv); %zscore

obj.ensembles.R = corr(deconv);