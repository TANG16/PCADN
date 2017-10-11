%% 
clear; close all;

load('/scratch/slim/klensink/data/pcadn/data.mat');

% Get the size
[nt nr ns] = size(D);

% Add Noise to the data
amp = 10;
skip = 2;
data_noise = addnoise(D,skip,amp);

% Take FFT
data_noise = fft(data_noise);
data_noise(end/2:end,:,:) = [];

% Spot operator to switch to midpoint offset
MH = opMH(nr,ns);
MO = nan(nt, nr, ns*2 -1);
MO_recon = nan(nt, nr, ns*2 -1);

% Number of modes to keep
k=801;

% Loop over all slices
for i=[20 400];
	
	% Take a time slice
	slice = squeeze(data_noise(i,:,:));
	
	% Switch to MO and store it
	MO(i,:,:) = reshape(MH*slice(:),401,2*401-1);

	% Perform PCA on slice and keep k modes
	
	[MO_recon(i,:,:), var(i,:)] = eig_pca(squeeze(MO(i,:,:)), k);
	
end %for

figure(1)
for i = [20 400]; 
	plot(cumsum(var(i,end:-1:1)/sum(var(i,:)))); 
	hold on; 
end

legend('f = 20','f = 400')
title('Cum Var TF - High Noise (10,2)')
xlabel('PCA Mode')
ylabel('Normalized Cumulative Variance')





