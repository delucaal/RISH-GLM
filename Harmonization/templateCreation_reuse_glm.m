%{

This program is created for dMRI Harmonization applying the method mentioned in:

Cross-site harmonization of diffusion MRI data without matched training subjects

 View ORCID ProfileAlberto De Luca, Tine Swartenbroekx,  View ORCID ProfileHarro Seelaar,  View ORCID ProfileJohn van Swieten,  View ORCID ProfileSuheyla Cetin Karayumak,  View ORCID ProfileYogesh Rathi,  View ORCID ProfileOfer Pasternak,  View ORCID ProfileLize Jiskoot,  View ORCID ProfileAlexander Leemans

https://www.biorxiv.org/content/10.1101/2024.05.01.591994v2

As the INPUTS, the user needs to set the following

SITES{1}.fullpath:           full path of the first site file in xls type
SITES{1}.name:               name of the first site
SITES{2}.fullpath:           full path of the second site file in xls type
SITES{2}.name:               name of the first site
templatesPath:               full path of the rish templates

OUTPUT:
HarmonizedDWIPaths:  the harmonized DWI subjects.

Author: Alberto De Luca, a.deluca-2@umcutrecht.nl
%}



function SITES=templateCreation_reuse_glm(SITES, option, templatesPath, existing_template )


mkdir(templatesPath);

%% USER INPUT -->
% Please enter the full path of xls files.
% Each file is prepared separetely for each site/scanner.
% File contains the fullpath of the dwi file and full path of brain mask in nhdr or nifti/analyze type

% images for template creation

%% Steps for harmonization:
% STEP1: Generate a case list of the images to be used
% STEP2: Save Ls (L0 L2 ...) of the above image.
% STEP3: Create rish feature templates

fprintf('STEP 1: Generating a case list of the images to be used\n');
SITES = loadData(SITES, option);


computeDTIMeasures(SITES, option);

fprintf('STEP 2: Saving Ls (L0 L2 ...) for each image\n');
SITES = computeRishFeatureImages(SITES, option);


if option.resample || option.bMap  || option.denoise
    computeDTIMeasures(SITES, option);
end


fprintf('STEP 3: Creating rish feature templates\n');
createTemplateBasedOnExisting_aluca_GLM(SITES,templatesPath,existing_template,option);
disp('End of template creation...');

% Clean the FSL folder
if ~option.debug
    cleanOutput(SITES);
end



end
