%% 
% PCA Denoise in Time Domain
clear; close all;

load('/scratch/slim/klensink/data/pcadn/data.mat');

% Get the size
[nt nr ns] = size(D);

% Add Noise to the data
amp = 10;
skip = 1;
data_noise = addnoise(D,skip,amp);

% Spot operator to switch to midpoint offset
MH = opMH(nr,ns);
MO = nan(nt, nr, ns*2 -1);
data_recon = zeros(size(D));
% Number of modes to keep
K=100;
modes = 14;
count = 1;


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
		
		% Store as data
		data_recon(i,:,:) = slice_recon;
		
		% Calc Errors
		RMSE(k,count) = sqrt(mean(mean((squeeze(D(i,:,:)) - 					slice_recon).^2)));
		SNR(k,count)  = -20*log10(norm(squeeze(D(i,:,:))-slice_recon,'fro')/norm(squeeze(D(i,:,:)),'fro'));
		
		
	end %for k

%imagesfigure(1)	
%plot(RMSE(:,count),'Color',[(nt-i)/nt, 0, 1-(nt-i)/nt]);
%hold on;

%figure(2)
%plot(SNR(:,count),'Color',[(nt-i)/nt, 0, 1-(nt-i)/nt]);
%hold on;
	
count = count+1;	
end %for i

%Slice to view
j = 10;

figure;
subplot 221
	imagesc(squeeze(D(:,:,j)))
	title('Original')
	colormap('gray')
	xlabel('Receiver')
	ylabel('Sample')
	caxis([-7 7]*1e1)
	
	set(gca,'Fontsize', 14, 'FontName', 'helvetica', 'FontWeight', 'demi');
subplot 222
	imagesc(squeeze(data_noise(:,:,j)))
	title('Added Noise')
	colormap('gray')
	xlabel('Receiver')
	ylabel('Sample')
	caxis([-7 7]*1e1)
	set(gca,'Fontsize', 14, 'FontName', 'helvetica', 'FontWeight', 'demi');
	
subplot 223
	imagesc(squeeze(data_recon(:,:,j)))	
	title('Recovered')
	colormap('gray')
	xlabel('Receiver')
	ylabel('Sample')
	caxis([-7 7]*1e1)
	set(gca,'Fontsize', 14, 'FontName', 'helvetica', 'FontWeight', 'demi');
	
subplot 224
	residual = D - data_recon;
	imagesc(residual(:,:,j))	
	title('Residual')
	colormap('gray')
	xlabel('Receiver')
	ylabel('Sample')
	caxis([-7 7]*1e1)
	set(gca,'Fontsize', 14, 'FontName', 'helvetica', 'FontWeight', 'demi');
		
	
	
	
	
	


