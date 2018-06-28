function [] = main()

if ~isdeployed
    disp('loading paths (HPC)')
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/jsonlab'))

    %addpath(genpath('/usr/local/jsonlab'))
    %addpath(genpath('/usr/local/vistasoft'))
end

% load my own config.json
config = loadjson('config.json');

normandflip(config)

dwi = niftiRead(config.dwi);

res = dwi.pixdim(1:3);
clear dwi

dwParams = dtiInitParams;
dwParams.clobber           =  1;
dwParams.eddyCorrect       = -1;
dwParams.phaseEncodeDir    = 2; 
dwParams.rotateBvecsWithRx = 0;
dwParams.rotateBvecsWithCanXform = 0;
dwParams.bvalsFile  = 'dwi.bvals';
% dwParams.outDir     = pwd;
dwParams.dwOutMm    = res;

%apply config params
if isfield(config, 'eddyCorrect')
    dwParams.eddyCorrect = str2double(config.eddyCorrect);
    dwParams.rotateBvecsWithRx = config.rotateBvecsWithRx;
    dwParams.rotateBvecsWithCanXform = config.rotateBvecsWithCanXform;
    dwParams.phaseEncodeDir    = str2num(config.phaseEncodeDir); 
end

dirs = {'noflip', 'xflip', 'yflip', 'zflip'};
bvecs_files = {'dwi_noflip.bvecs', 'dwi_xflip.bvecs', 'dwi_yflip.bvecs',...
        'dwi_zflip.bvecs'};
trilin_name = {'noflip_dti_trilin', 'xflip_dti_trilin', 'yflip_dti_trilin',...
        'zflip_dti_trilin'};
fg_names = {'noflip_fg', 'xflip_fg', 'yflip_fg', 'zflip_fg'};

    
for i = 1:length(dirs)
	mkdir(char(dirs(i)))
end    
long_fibers = zeros(1,length(dirs));

%% params

opts.faMaskThresh = 0.30;

% Distance between steps in the tractography algoithm
opts.stepSizeMm = 1;
% Stopping criteria FA<0.2
opts.faThresh = 0.2;
% Discard Fibers shorter than 50mm or longer than 250mm
opts.lengthThreshMm = [50 250];
% Stopping criteria angle between steps >30 degrees
opts.angleThresh = 30;
% Unknown.....
opts.wPuncture = 0.2;
% There are multip.e algorithms that can be used.  We prefer STT. See:
% Basser PJ, Pajevic S, Pierpaoli C, Duda J, Aldroubi A. 2000.
% In vivo fiber tractography using DT-MRI data.
% Magnetic Resonance in Medicine 44(4):625-32.
opts.whichAlgorithm = 1;
% Interpolation method. After each step we interpolate the tensor at that
% point. Trilinear interpolation works well.
opts.whichInterp = 1;
% This adds some randomness to each seed point. Each seed point is move
% randomly by randn*.1mm
opts.offsetJitter = 0;
% We seed in voxel in multiple locations. [0.25 and 0.75] Seeds each voxel
% at 8 equidistant locations.  For one seed in the middle of the voxel use

opts.seedVoxelOffsets = [0.25 0.75];


%% loop through flips
for i = 1:length(bvecs_files)
    dwParams.outDir = strcat('./', char(dirs(i)));
    dwParams.dt6BaseName = char(trilin_name(i));
    dwParams.bvecsFile  = char(bvecs_files(i));
    dtiInit(config.dwi, config.t1, dwParams)
    
    [dt6, xformToAcpc, mmPerVox] = dtiLoadTensorsFromNifti(...
        strcat('./',char(dirs(i)), '/', char(trilin_name(i)),'/bin/tensors.nii.gz'));
    % Compute FA at every voxel
    fa = dtiComputeFA(dt6);
    % Sometimes noise can cause impossible FA values so we clip them
    fa(fa>1) = 1; fa(fa<0) = 0;
    
    
    %% Create an ROI for wholebrain tractography
    roiAll = dtiNewRoi('all');
    % Make a mask of voxels with FA>faMaskThresh
    mask = fa >= opts.faMaskThresh;
    % Convert mask image to a list of coordinates
    [x,y,z] = ind2sub(size(mask), find(mask));
    clear mask fa;
    % Transofrm mask coordinates to subjects ACPC space
    roiAll.coords = mrAnatXformCoords(xformToAcpc, [x,y,z]);
    clear x y z;
    % Smooth the ROI and fill holes
    roiAll = dtiRoiClean(roiAll,3,{'fillHoles'});
    
    %% Perform wholebrain tractography
    fg = dtiFiberTrack(dt6, roiAll.coords, mmPerVox, xformToAcpc, 'wholeBrain', opts);
    save(strcat(char(fg_names(i)), '.mat'),'fg','-v7.3');
    %dtiExportFibersMrtrix(fg, strcat(char(fg_names(i)), '.tck'));
    
    numlongfibers = 0;
    for istreamlines=1:length(fg.fibers)
        if length(fg.fibers{istreamlines}) > 50
            numlongfibers = numlongfibers + 1;
        end
    end
    long_fibers(i) = numlongfibers;
end

[M, I] = max(long_fibers);

flipdirection = {'no flip', 'x flip', 'y flip' , 'z flip'};

fileID = fopen('results.txt', 'w');
fprintf(fileID, 'Gradient Flip Recommended: %s \n', char(flipdirection(I)));
fprintf(fileID, 'Number of long fibers:\n');
for line = 1:length(long_fibers)
    fprintf(fileID, '%s: %d \n', char(flipdirection(line)), long_fibers(line));
end
% fprintf(fileID, 'no flip: %d \n', long_fibers(1));
% fprintf(fileID, 'x flip: %d \n', long_fibers(2));
% fprintf(fileID, 'y flip: %d \n', long_fibers(3));
% fprintf(fileID, 'z flip: %d \n', long_fibers(4));
fclose(fileID);

%% product.json generation

message_item = struct;
message_item.type = 'info';
% 32 used to represent a space (since ' ' doesn't work)
message_item.msg = strcat(flipdirection{I}, 32, 'recommended');

plot_item = struct;
plot_item.type = 'plotly';
plot_item.name = 'Flip Recommendation';
plot_item.data = struct;
plot_item.layout = struct;

plot_item.data.values = long_fibers;
plot_item.data.labels = flipdirection;
plot_item.data.type = 'pie';
plot_item.data = { plot_item.data };

results = { message_item, plot_item };
savejson('brainlife', results, 'product.json');

end

