function HarmonizedDWIPaths=harmonization(xlsFilePath_REF, REF_name, xlsFilePath_OTHER, OTHER_name, option, templatesPath )

SITES{1}.fullpath = xlsFilePath_REF;  SITES{1}.name=REF_name; %REFERENCE SITE
SITES{2}.fullpath = xlsFilePath_OTHER;  SITES{2}.name=OTHER_name; % TARGET SITE
option.name = OTHER_name;
%% Steps for harmonization:

%% STEP1: Load data to be harmonized
fprintf('STEP 1: Generating a case list of the images to be used\n');
SITES = loadData(SITES, option);



%% STEP 2: if option.FA = 1, compute dti measures for the original data
fprintf('STEP 2: Computing DTI measures\n');
computeDTIMeasures(SITES, option);


%% STEP 3: compute RISH features
fprintf('STEP 3: Saving Ls (L0 L2 ...) for each image\n');
SITES = computeRishFeatureImages(SITES, option);


if option.resample || option.bMap  || option.denoise
    computeDTIMeasures(SITES, option);
end

%% Optional steps: debugging mode
if option.debug
    fprintf('Optional: Registering original data to MNI Space...\n');
    WarpOriginalDiffusionMeasuresToMNISpace(SITES, option);
    
    fprintf('Optional: Registering resampled data to template space...\n');
    WarpDiffusionMeasuresToHarmonizationTemplate(SITES, templatesPath, option);
    
    fprintf('Optional: Registering the data in the template space to MNI space...\n');
    WarpHarmonizationTemplateToMNISpace(SITES, templatesPath, option);
end

%% STEP 4: Harmonization
fprintf('STEP 4: Harmonizing data\n');
[HarmonizedDWIPaths, SITES] = harmonize_CleanOutliers(SITES,templatesPath, option);

option.siteno=2; % apply the rest only to the target site
%% Optional steps: debugging mode
if option.debug
    fprintf('Optional: Computing DTI measures for harmonized data...\n');
    computeDTIMeasures(SITES, option);
    
    fprintf('Optional: Warp harmonized data to template Space...\n');
    WarpDiffusionMeasuresToHarmonizationTemplate(SITES, templatesPath, option);
    
    fprintf('Optional: Warp harmonized data to MNI Space...');
    WarpHarmonizationTemplateToMNISpace(SITES, templatesPath, option);
end

% Clean the FSL folder
if ~option.debug
    cleanOutput(SITES);
end
end