clc, clear
addpath(genpath('./Harmonization'));

%% REQUIRED INPUT:
%Harmonization reference site, templateCreation input:  [filePath_REF, '_dwi.txt'] and
%[filePath_REF, '_mask.txt']  keep the paths of the dwi paths
%and the mask paths, but only prefix is given as input. See the example
%below.
SITES{1}.fullpath = 'Site1'; % Will try to read files named 'Site1_dwi.txt' and 'Site1_mask.txt' containing a list of diffusion data and corresponding masks
SITES{1}.name = 'Site1'; % A short acronym for this site

SITES{2}.fullpath = 'Site2'; % same as for site1
SITES{2}.name = 'Site'; % same as for site1

SITES{3}.fullpath = 'Site3'; % same as for site1
SITES{3}.name = 'Site3'; % same as for site1

%path of the template  to be created
templatesPath = 'Path to a RISH-GLM template. It will be created by this script if it does not exist.';
baseTemplate = 'Path to a RISH template. It will be created if it does not exist (equivalent to example_main_rish.m)';

%Creates the option.harmonizeDir under each subject directory, stores the
%harmonization results inside this directory.
option.covariates = load('GLM_Covariates.txt');

option.name='AIO';

option.robust_std = 0; % No robust average
option.reference_site = 2; % Which site to use as a reference
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

option.harmonizeDir = 'harmonizedOutputRISH'; 
SITES=templateCreation(filePath_REF, REF_name, filePath_OTHER, OTHER_name, option, templatesPath);
option.harmonizeDir = 'harmonizedOutputRISHGLM'; 
SITES = templateCreation_reuse_glm(SITES, option, templatesPath, baseTemplate);
% return
%% PART 2: Harmonization: Template creation should be succesfully run to move to the harmonization part

SITES{1}.fullpath = 'Source_LEID_1';
SITES{1}.name = 'R4_8ch';

SITES{2}.fullpath = 'Source_LEID_3_4';
SITES{2}.name = 'R5_8_32ch';

SITES{3}.fullpath = 'Source_RTM_5_6';
SITES{3}.name = 'RTM_48ch';

option.covariates = load('GLM_Covariates.txt');

HarmonizedDWIPaths=harmonization_touchRef_glm(SITES, option, templatesPath );
