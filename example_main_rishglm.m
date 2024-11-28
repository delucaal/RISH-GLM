% clc, clear
% addpath(genpath('./Utilities'));
% addpath(genpath('./NIfTI_20140122'));

%% REQUIRED INPUT:
%Harmonization reference site, templateCreation input:  [filePath_REF, '_dwi.txt'] and
%[filePath_REF, '_mask.txt']  keep the paths of the dwi paths
%and the mask paths, but only prefix is given as input. See the example
%below.
SITES{1}.fullpath = 'Source_LEID_1';
SITES{1}.name = 'R4_8ch';

SITES{2}.fullpath = 'Source_LEID_3_4';
SITES{2}.name = 'R5_8_32ch';

SITES{3}.fullpath = 'Source_RTM_5_6';
SITES{3}.name = 'RTM_48ch';

%path of the template  to be created
templatesPath='/autofs/arch11/DATA/PROVIDI_LAB/Alexander/My_data/FTD_RISC/ProcessedData/Harmonization/Training_GLM_AIO/TemplateHarmGLM_RefS2Z_NoS3_Fix2/';
baseTemplate = '/autofs/arch11/DATA/PROVIDI_LAB/Alexander/My_data/FTD_RISC/ProcessedData/Harmonization/Training_GLM_AIO/TemplateHarmRegistrationTest/';

%Creates the option.harmonizeDir under each subject directory, stores the
%harmonization results inside this directory.
option.harmonizeDir = 'harmonizedOutputGLM_WCovsRefS2ZNoS3_Fix'; 
option.covariates = load('GLM_Covariates_Fix.txt');
option.covariates = option.covariates;%.^(2);
for x=1:size(option.covariates,2)
   if(max(option.covariates(:,x)) > 1)
       V = option.covariates(:,x);
       V = (V-mean(V))/std(V);
       option.covariates(:,x) = V;
   end
end

option.name='AIO';

option.robust_std = 0; % No robust average
option.reference_site = 2; % Site 2 is the reference
option.smooth_std = 0; % no smoothing

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
SITES = templateCreation_reuse_glm(SITES, option, templatesPath, baseTemplate);
% return
%% PART 2: Harmonization: Template creation should be succesfully run to move to the harmonization part
% filePath_REFharmonization='G1'; REF_name='R5_8_32ch'; % 8ch and 32ch are inverted... 
% filePath_OTHERharmonization='G2'; OTHER_name='R4_8ch'; %  list of all subjects in target site to be harmonized 

SITES{1}.fullpath = 'Source_LEID_1';
SITES{1}.name = 'R4_8ch';

SITES{2}.fullpath = 'Source_LEID_3_4';
SITES{2}.name = 'R5_8_32ch';

SITES{3}.fullpath = 'Source_RTM_5_6';
SITES{3}.name = 'RTM_48ch';

load('/autofs/arch11/DATA/PROVIDI_LAB/Alexander/My_data/FTD_RISC/ProcessedData/DEMOGRAPHICS_MATCHED_INC_PRECONSYM_NODUP_WEC_AGEC2_TS.mat');

option.covariates = cell(length(SITES),1);
option.covariates_2_correct = []; %NOF SITES + the actual index. For now not used

for s=1:length(SITES)
    list = read_files_list([SITES{s}.fullpath '_dwi.txt']);
    option.covariates(s) = {[]};
    
    for ix=1:length(list)
        [fp,fn,ext] = fileparts(list{ix});
        for fu=0:6
            fn = strrep(fn,['FU' num2str(fu)],['FU_' num2str(fu)]);
        end

        idx = locateindb(DEMOGRAPHICS_MATCHED,fn);
        sex = strtrim(DEMOGRAPHICS_MATCHED{idx,6});

        if(strcmpi(sex,'man') || strcmpi(sex,'male'))
            sex = 1;
        elseif(strcmpi(sex,'vrouw') || strcmpi(sex,'female'))
            sex = 0;
        else
            error('unexpected');
        end
        age = DEMOGRAPHICS_MATCHED{idx,8};
        if(isempty(age))
            age = 0;
        end
        site = DEMOGRAPHICS_MATCHED{idx,17};
        option.covariates{s} = cat(1,option.covariates{s},[age sex]); % other coil + age + sex
    end

end

HarmonizedDWIPaths=harmonization_touchRef_glm(SITES, option, templatesPath );

%% Auxiliary function for covariates
function idx = locateindb(db,name)
    idx = -1;
    stp = strfind(name,'_dwi_FP');
    std = strfind(name,'sub-');
    name = name(std:stp-1);
    for ix=2:size(db)
       if(strcmp(db{ix,2},name))
           idx = ix;
           return;
       end
    end
end

function list = read_files_list(file)
    f = fopen(file,'rt');
    list = {};
    while(~feof(f))
        l = fgetl(f);
        list = cat(1,list,{l});
    end
    fclose(f);
end
