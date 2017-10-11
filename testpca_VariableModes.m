function [slice_recon, var] = testpca(slice, k)
%% Take PCA of a 2D image
% Returns the image reconstructed out of the first k modes found
% using eig PCA, and the variance explained by each mode
% 
% Use: [slice_recon, var] = testpca(slice, k);


% Get slice Dimensions
[nr ns] = size(slice);

% De-Mean the columns of the image
col_mean = repmat(mean(slice, 1),nr,1);
slice_dm = slice - col_mean;

% Get cov matrix
C = cov(slice_dm);

% Get Eigenvalues and eigenvectors
[eigenvectors eigenvalues] = eig(C);

% eigenvectors has each eigenvector in a column
% Keeping 'k' eigenvectors
eigenvectors = eigenvectors(:,end-k+1:end);

slice_final = eigenvectors' * slice_dm';

slice_recon = (eigenvectors*slice_final)' + col_mean;


