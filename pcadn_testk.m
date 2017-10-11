%% 
clear; close all;

load('/scratch/slim/klensink/data/pcadn/data.mat');

% Get the size
[nt nr ns] = size(D);

% Add Noise to the data
amp = 10;
skip = 10;
data_noise = addnoise(D,skip,amp);

% Spot operator to switch to midpoint offset
MH = opMH(nr,ns);
MO = nan(nt, nr, ns*2 -1);
MO_recon = nan(nt, nr, ns*2 -1);

% Number of modes to keep
K=100;
modes = 52;
count = 1;

figure(1)

% Loop over slices
for i=1:nt;

	% Loop over all modes
	for k = modes;
	
		% Take a time slice
		slice = squeeze(data_noise(i,:,:));
	
		% Switch to MO
		MO = reshape(MH*slice(:),401,2*401-1);

		% Perform PCA on slice and keep k modes	
		[MO_recon, var(i,:)] = eig_pca(MO, k);
	
		% Take the adjoint opMH and return to aq domain
		slice_recon = reshape(MH'*MO_recon(:),401,401);
	
		RMSE(k,count) = sqrt(mean(mean((squeeze(D(i,:,:)) - 					slice_recon).^2)));
		SNR(k,count)  = -20*log10(norm(squeeze(D(i,:,:))-slice_recon,'fro')/norm(squeeze(D(i,:,:)),'fro'));
		
		
		
	end %for k

figure(1)	
plot(RMSE(:,count),'Color',[(nt-i)/nt, 0, 1-(nt-i)/nt]);
hold on;

figure(2)
plot(SNR(:,count),'Color',[(nt-i)/nt, 0, 1-(nt-i)/nt]);
hold on;
	
count = count+1;	
end %for i

title('RMSE of PCA applied to time slices. Depth ranges from shallow (red) to deep (blue)')
xlabel('Number of PCA Modes')
ylabel('RMSE')

[y i] = min(RMSE);

min_error = mean(y)

figure(2)
plot(i)


