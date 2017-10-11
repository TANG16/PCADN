function [data] = addnoise(data, skip, amplitude);
%% Add noise to every skipth row of data
% Use: [data] = addnoise(data, skip, amplitude)

[nt nr ns] = size(data);

ind = 1:skip:ns;

% Add Noise
data(:,ind,:) = data(:,ind,:) + amplitude*randn(nt, length(ind), ns);
