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



function SITES=templateCreation(SITES, option, templatesPath )


mkdir(templatesPath);

% Addition: automatically take care of logistics - A. De Luca (a.deluca-2@umcutrecht.nl)
for s=1:length(SITES)
    if(exist([SITES{s}.fullpath filesep 'desc_dwi.txt'],'file') < 1)
        dwis = dir(fullfile(SITES{s}.fullpath,'*.bval'));
        dwis_names = cell(length(dwis),1);
        mask_names = cell(length(dwis),1);

        for d=1:length(dwis)
            bn = strrep(dwis(d).name,'.bval','');
            dn = dir(fullfile(dwis(d).folder,[bn '.nii*']));
            mn = dir(fullfile(dwis(d).folder,[bn '*mask*.nii*']));
            if(length(dn) ~= 1 || length(mn) ~= 1)
                error(['Either no or multiple hits for ' bn]);
            end
            dwis_names{d} = fullfile(dn.folder,dn.name);
            mask_names{d} = fullfile(mn.folder,mn.name);
        end

        fp = fopen([SITES{s}.fullpath filesep 'desc_dwi.txt'],'wt');
        fm = fopen([SITES{s}.fullpath filesep 'desc_mask.txt'],'wt');
        for d=1:length(dwis_names)
            fprintf(fp,'%s%s',dwis_names{d},newline);
            fprintf(fm,'%s%s',mask_names{d},newline);
        end
        fclose(fp);
        fclose(fm);
    end
    SITES{s}.fullpath = [SITES{s}.fullpath filesep 'desc'];
end

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
createTemplates(SITES,templatesPath,option);
disp('End of template creation...');

% Clean the FSL folder
if ~option.debug
    cleanOutput(SITES);
end



end
