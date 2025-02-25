clc, clear
addpath(genpath('./Harmonization'));
mkdir('log');
%% REQUIRED INPUT:
%Harmonization reference site, templateCreation input:  [filePath_REF, '_dwi.txt'] and
%[filePath_REF, '_mask.txt']  keep the paths of the dwi paths
%and the mask paths, but only prefix is given as input. See the example
%below.
SITES{1}.fullpath = 'to_be_completed'; % A folder containing the training data of Site 1: diffusion MRI data (.nii), a mask (same name as data + _mask.nii), and the corresponding .bval/.bvec. 
SITES{1}.name = 'name_of_site1';

SITES{2}.fullpath = 'to_be_completed'; % A folder containing the training data of Site 2: diffusion MRI data (.nii), a mask (same name as data + _mask.nii), and the corresponding .bval/.bvec. 
SITES{2}.name = 'name_of_site2';

% Note: you need at least 2 sites to run the harmonization.

SITES{3}.fullpath = 'to_be_completed'; % A folder containing the training data of Site 3: diffusion MRI data (.nii), a mask (same name as data + _mask.nii), and the corresponding .bval/.bvec. 
SITES{3}.name = 'name_of_site3';

% Note: More sites can be added. 

%path of the template  to be created
templatesPath = 'to_be_completed';
baseTemplate = 'to_be_completed';

% Covariates should be specified in a text file using one column per covariate of interest.
option.covariates = load('path_to_covariates.txt');

option.name='AIO';

option.robust_std = 0; % No robust average. Can be useful if the data contains outliers.
option.reference_site = 1; % Which site to use as a reference. By default, this is the first site.
option.smooth_std = 0; % no smoothing

% Perform Z-score normalization of all covariates
for x=1:size(option.covariates,2)
   if(max(option.covariates(:,x)) > 1)
       V = option.covariates(:,x);
       V = (V-mean(V))/std(V);
       option.covariates(:,x) = V;
   end
end

%% OTHER OPTIONS:
option.force = 0; % forces overwriting existing files when option.force is 1, default is 0
option.denoise = 0; % set to 1 to denoise dwi data, default is 0
option.debug = 0; % when option.debug = 1, it will store all  diffusion measure files in subject space, MNI space and template space.
% if you don't want to run this optopn, set option.debug = 0, it is going to remove the FSL folder at the end of harmonization, default is 1.

% set SHOrder=2 at least for 6 gradient directions, SHOrder=4 at least for 15
% gradient directions, SHOrder=6  at least for 28 gradient directions,
% SHOrder=8 for at least 45 gradient directions, default is 6.
option.SHOrder = 6;

option.bMap = 0; % map b val when value is 1, default is 0
option.bValue =1000; % new b value: this option works only when option.bMap=1

option.resample = 0; % resample when value is 1, default is 0
option.sp_high = [1.7188,1.7187,3];% new spatial resolution: this option works only when option.resample=1

option.travelingHeads = 0; %set it to 1  when harmonizing traveling heads (the same subjects scanned using different scanners). 

option.harmonizedName = 'harmonized_'; % prefix for harmonized dwi data output
option.rishedName = 'rished_'; % prefix for non-harmonized dwi data output. This data has been interpolated using spherical harmonics and provides a fair comparison to evaluate harmonization effects.
option.DTIdir='/DTI/'; % diffusion measures output folder name, it is created as a subfolder under each subject directory

% option.lambda=0.006; % regularization parameter for SH. To be adjusted
% per dataset. Higher -> more regularization = more loss of data
% specificity. Lower -> more sensitivity to noise and unstable
% harmonization
option.lambda = 0.001;
%% START OF THE PROGRAM:

%% PART 1: Creating Harmonization Templates
option.parallel=0; %parallel computing for ants  
option.numcores=4;% number of cpu cores
option.siteno=1;% create RISH features and diffusion measures for each subject, don't change this!

option.harmonizeDir = 'harmonizedOutputRISH'; % creates a folder within the source of each site to store the harmonized data (conventional RISH)
SITES = templateCreation(SITES, option, baseTemplate);
option.harmonizeDir = 'harmonizedOutputRISHGLM'; % creates a folder within the source of each site to store the harmonized data (RISH-GLM)
SITES = templateCreation_reuse_glm(SITES, option, templatesPath, baseTemplate);

%% PART 2: Harmonization: Template creation should be succesfully run to move to the harmonization part

% PART 1 was about training harmonization. 
% This part follows the same structure but allows to apply harmonization to new data.
SITES{1}.fullpath = 'to_be_completed';
SITES{1}.name = 'name_of_site1';

SITES{2}.fullpath = 'to_be_completed';
SITES{2}.name = 'name_of_site2';

SITES{3}.fullpath = 'to_be_completed';
SITES{3}.name = 'name_of_site3';

option.covariates = load('to_be_completed.txt');
option.covariates_2_correct = []; %NOF SITES + the actual index. Not used for now
HarmonizedDWIPaths=harmonization_touchRef_glm(SITES, option, templatesPath );
