%% 
% PCA Denoise in Frequency Domain
clear; close all;

load('/scratch/slim/klensink/data/pcadn/data.mat');

% Get the size
[nt nr ns] = size(D);

tic;

% Add Noise to the data
amp = 10;
skip = 1;
data_noise = addnoise(D,skip,amp);

% Take 1D FFT in time SKIP NOISE RN
data_f = fft(data_noise,[],1);
keep = data_f;

% Keep only non - negative frequencies
data_f = data_f(1:end/2+1,:,:); 
[nf nr ns] = size(data_f);

% Setup of spot operator to switch to midpoint offset
MH = opMH(nr,ns);

% Initialize empty frequency volume for denoising results before ifft
% as well as volume for reconstructed data
f_vol = zeros(nt, nr, ns);
data_recon = zeros(size(D));
count = 1;

% Number of modes to keep 
modes = 1:100;

% Loop over all modes
for k = modes;
	
	% Loop over slices
	for i=1:nf;
	
		% Take a frequency slice
		slice = squeeze(data_f(i,:,:));
	
		% Switch to MO
		MO = reshape(MH*slice(:),401,2*401-1);

		% Perform PCA on slice and keep k modes	
		[MO_recon, var(i,:)] = eig_pca(MO, k);
	
		% Take the adjoint opMH and return to aq domain
		slice_recon = reshape(MH'*MO_recon(:),401,401);
		
		% Store into denoised frequency slice volume
		f_vol(i,:,:) = slice_recon;			
					
		
	end %for i


	% Take complex conj and reconstruct negative frequencies
	f_vol(nf+1:end,:,:) = conj(flipud(f_vol(2:nf-1,:,:)));

	% 1D IFFT to get back reconstructed data
	data_recon = ifft(f_vol,[],1);



	% Calulate RMSE
	for ii = 1:ns
		RMSE(ii,count) = sqrt(mean(mean((D(:,:,ii) - data_recon(:,:,ii)).^2)));
	end % For ns
	
	count = count+1;
end %for k

%Slice to view
j = 10;
	
figure;
subplot 131
	imagesc(squeeze(D(:,:,j)))
	title('Original')
	colormap('gray')
	xlabel('Receiver')
	ylabel('Sample')
	caxis([-7 7]*1e1)
	
	set(gca,'Fontsize', 14, 'FontName', 'helvetica', 'FontWeight', 'demi');
subplot 132
	imagesc(squeeze(data_noise(:,:,j)))
	title('Added Noise')
	colormap('gray')
	xlabel('Receiver')
	ylabel('Sample')
	caxis([-7 7]*1e1)
	set(gca,'Fontsize', 14, 'FontName', 'helvetica', 'FontWeight', 'demi');
	
subplot 133
	imagesc(squeeze(data_recon(:,:,j)))	
	title('Recovered')
	colormap('gray')
	xlabel('Receiver')
	ylabel('Sample')
	caxis([-7 7]*1e1)
	set(gca,'Fontsize', 14, 'FontName', 'helvetica', 'FontWeight', 'demi');
	
figure;
plot(mean(RMSE))	

save('/home/slim/klensink/figures/test/pcadn/RMSE_1_100.mat','RMSE','modes')

toc	
	
	
	
	


