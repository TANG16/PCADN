clear; close all;

load('/scratch/slim/klensink/data/slice_tx.mat');



% De-Mean the columns of the image
col_mean = repmat(mean(slice, 1),401,1);
slice_dm = slice - col_mean;

% Get cov matrix
C = cov(slice_dm);

% Get Eigenvalues and eigenvectors
[eigenvectors eigenvalues] = eig(C);

% eigenvectors has each eigenvector in a column
% Keeping 'k' eigenvectors
k = 30;
eigenvectors = eigenvectors(:,end-k+1:end);

slice_final = eigenvectors' * slice_dm';

slice_recon = (eigenvectors*slice_final)' + col_mean;

% Original Image
figure(1)
subplot 211
imagesc(slice);
title('Original Image')
colormap('gray')

subplot 212
imagesc(slice_recon);
title('Reconstructed Image')
colormap('gray')

