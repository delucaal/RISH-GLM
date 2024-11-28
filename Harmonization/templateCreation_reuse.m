%{

This program is created for dMRI Harmonization applying the mothod mentioned in:

Multi site Harmonization of diffusion MRI data in a registration framework (ANTs based).

If you use this code, please cite the following:
Hengameh Mirzaalian, Lipeng Ning, Peter Savadjiev, Ofer Pasternak, Sylvain Bouix, Oleg Michailovich,
Gerald Grant, Christine E. Marx, Rajendra A. Morey, Laura A. Flashman, Mark S. George, Thomas. McAllister,
Norberto Andaluz, Lori Shutter, Raul Coimbra, Ross D. Zafonte11, Mike J. Coleman, Marek Kubicki,
Carl-Fredrik Westin, Murray B. Stein, Martha E. Shenton, and Yogesh Rathi. 2016.
 “Multi-site harmonization of diffusion MRI data in a registration framework.”
Brain Imaging and Behavior Journal.

As the INPUTS, the user needs to set the following

SITES{1}.fullpath:           full path of the first site file in xls type
SITES{1}.name:               name of the first site
SITES{2}.fullpath:           full path of the second site file in xls type
SITES{2}.name:               name of the first site
templatesPath:               full path of the rish templates

e.g:
SITES{1}.fullpath = 'site1.xls'; SITES{1}.name='SITE1';
SITES{2}.fullpath = 'site2.xls'; SITES{2}.name='SITE2';
templatesPath='/projects/schiz/suheyla/suheyla-harmonization/Templates1/'; mkdir(templatesPath);

OUTPUT:
HarmonizedDWIPaths:  the harmonized DWI subjects.

Author: Suheyla Cetin Karayumak, skarayumak@bwh.harvard.edu
%}



function SITES=templateCreation_reuse(filePath_REF, REF_name, filePath_OTHER, OTHER_name, option, templatesPath, existing_template )


mkdir(templatesPath);

%% USER INPUT -->
% Please enter the full path of xls files.
% Each file is prepared separetely for each site/scanner.
% File contains the fullpath of the dwi file and full path of brain mask in nhdr or nifti/analyze type

% images for template creation
SITES{1}.fullpath = filePath_REF;  SITES{1}.name=REF_name; %REFERENCE SITE
SITES{2}.fullpath = filePath_OTHER;  SITES{2}.name=OTHER_name;
option.name=OTHER_name;

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
createTemplateBasedOnExisting_aluca(SITES,templatesPath,existing_template,option);
disp('End of template creation...');

% Clean the FSL folder
if ~option.debug
    cleanOutput(SITES);
end



end
