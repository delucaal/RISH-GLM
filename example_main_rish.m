% clc, clear
% addpath(genpath('./Utilities'));
% addpath(genpath('./NIfTI_20140122'));

%% REQUIRED INPUT:
%Harmonization reference site, templateCreation input:  [filePath_REF, '_dwi.txt'] and
%[filePath_REF, '_mask.txt']  keep the paths of the dwi paths
%and the mask paths, but only prefix is given as input. See the example
%below.
filePath_REF='Source_LEID_3_4'; REF_name='R5_8_32ch';
%Harmonization target site   templateCreation input.
filePath_OTHER='Source_RTM_5_6'; OTHER_name='RTM_48ch';

%path of the template  to be created
templatesPath='/autofs/arch11/DATA/PROVIDI_LAB/Alexander/My_data/FTD_RISC/ProcessedData/Harmonization/Training_GLM_AIO/TemplateHarm_CheckScales3456/';
baseTemplate = '/autofs/arch11/DATA/PROVIDI_LAB/Alexander/My_data/FTD_RISC/ProcessedData/Harmonization/Training_GLM_AIO/TemplateHarmRegistrationTest';

%Creates the option.harmonizeDir under each subject directory, stores the
%harmonization results inside this directory.
option.harmonizeDir = 'harmonizedOutputForCheckScales'; 


%% OTHER OPTIONS:
option.force = 0; % forces overwriting existing files when option.force is 1, default is 0
option.denoise = 0; % set to 1 to denoise dwi data, default is 0
option.debug = 1; % when option.debug = 1, it will store all  diffusion measure files in subject space, MNI space and template space.
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
option.rishedName = 'rished_'; % prefix for harmonized dwi data output
option.DTIdir='/DTI/'; % diffusion measures output folder name, it is created as a subfolder under each subject directory
% option.lambda=0.006; % regularization parameter for SH 
% Utrecht_Harmonization_Training
option.lambda = 0.001;
%% START OF THE PROGRAM:

%% PART 1: Creating Harmonization Templates
option.parallel=0; %parallel computing for ants  
option.numcores=4;% number of cpu cores
option.siteno=1;% create RISH features and diffusion measures for each subject, don't change this!
% SITES=templateCreation(filePath_REF, REF_name, filePath_OTHER, OTHER_name, option, templatesPath);
% SITES=templateCreation(filePath_REF, REF_name, filePath_OTHER, OTHER_name, option, templatesPath);
SITES = templateCreation_reuse(filePath_REF, REF_name, filePath_OTHER, OTHER_name, option, templatesPath, baseTemplate);
return
%% PART 2: Harmonization: Template creation should be succesfully run to move to the harmonization part
filePath_REFharmonization='./G1'; REF_name='R5_8_32ch'; 
% filePath_OTHERharmonization='./BELGIUM1_b0b1000_spatialtemplate'; OTHER_name='BELGIUM1'; %  list of all subjects in target site to be harmonized 
filePath_OTHERharmonization='./G2'; OTHER_name='R4_8ch'; %  list of all subjects in target site to be harmonized 
if exist([templatesPath  'Mean_' SITES{1}.name '_FA.nii.gz'], 'file') &&  exist([templatesPath  'Mean_' SITES{2}.name '_FA.nii.gz'], 'file')
	HarmonizedDWIPaths=harmonization_touchRef(filePath_REFharmonization, REF_name, filePath_OTHERharmonization, OTHER_name, option, templatesPath );
else
   error(['Template does not exist, check your templates folder first and be sure that ', [templatesPath  'Mean_' SITES{1}.name '_FA.nii.gz'], ' and ', ...
 [templatesPath  'Mean_' SITES{2}.name '_FA.nii.gz'], 'exist']);
end

