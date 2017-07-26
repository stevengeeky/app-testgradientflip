function [] = normandflip(config)

% normalizes the bvals and flips the bvecs


% Parameters used for normalization
params.thresholds.b0_normalize    = 200;
params.thresholds.bvals_normalize = 100;

%% Normalize HCP files to the VISTASOFT environment

bvals.val = dlmread(config.bvals);

% Round the numbers to the closest thousand 
% This is necessary because the VISTASOFT software does not handle the B0
% when they are not rounded.
[bvals.unique, ~, bvals.uindex] = unique(bvals.val);

bvals.unique(bvals.unique <= params.thresholds.b0_normalize) = 0;
bvals.unique  = round(bvals.unique./params.thresholds.bvals_normalize) ...
    *params.thresholds.bvals_normalize;
bvals.valnorm = bvals.unique( bvals.uindex );
dlmwrite('dwi.bvals',bvals.valnorm);

%% Flip the Bvecs on chosen axis

%load bvecs
bvecs = dlmread(config.bvecs);

if ~(size(bvecs,1) == 3), bvecs = bvecs'; end

dlmwrite('dwi_noflip.bvecs', bvecs);

bvecs_xflip = bvecs;
bvecs_xflip(1,:) = -bvecs(1,:);
dlmwrite('dwi_xflip.bvecs', bvecs_xflip);

bvecs_yflip = bvecs;
bvecs_yflip(2,:) = -bvecs(2,:);
dlmwrite('dwi_yflip.bvecs', bvecs_yflip);

bvecs_zflip = bvecs;
bvecs_zflip(3,:) = -bvecs(3,:);
dlmwrite('dwi_zflip.bvecs', bvecs_zflip);

end

