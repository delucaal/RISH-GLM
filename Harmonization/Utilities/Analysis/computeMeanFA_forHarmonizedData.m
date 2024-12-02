function SITES=computeMeanFA_forHarmonizedData(filePath_REF, REF_name, filePath_OTHER, OTHER_name, harmonizedDir, atlaspath)
setPaths;

SITES{1}.fullpath = filePath_REF;  SITES{1}.name=REF_name; %REFERENCE SITE
SITES{2}.fullpath = filePath_OTHER;  SITES{2}.name=OTHER_name;
option.harmonizeDir = harmonizedDir; %prisma st ge st
option.name=OTHER_name;

[SITES, err] = loadData(SITES, option);
atlas = load_untouch_nii(atlaspath);


s=1;
muAllorig=[];
for i=1:SITES{s}.ImageNum
    data  = SITES{s}.InputImages{i};
    
    load([data.Image_Harmonized_dir,'newfile_suffix.mat']);
    load([data.Image_Harmonized_dir,'newmask_suffix.mat']);
    load([data.Image_Harmonized_dir,'processes.mat']);
    SITES{s}.InputImages{i}.DWISuffix = newfile_suffix;
    SITES{s}.InputImages{i}.MaskSuffix = newmask_suffix;
    SITES{s}.InputImages{i}.denoised = process(1);
    SITES{s}.InputImages{i}.bValueMapped = process(2);
    SITES{s}.InputImages{i}.resampled = process(3);
    
    
    data  = SITES{s}.InputImages{i};
    % data.DWISuffix: can be '_rescaled', '_bmapped', or '_denoised'
    if ~isempty(data.DWISuffix)
        data.ImageName = [data.ImageName, data.DWISuffix];
        data.MaskName = [data.MaskName, data.MaskSuffix];
        data.MaskFullPath =  [data. Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
        data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
    end
    
    dataImg = [data.ImageDirectory 'FSL/' data.ImageName  '_InMNI_FA.nii.gz'];
    
    niiOrig = load_untouch_nii(dataImg);
    
    
    ind=atlas.img == 1;
    wholeOrigWm=niiOrig.img.*single(atlas.img);
    
    
    muAllorig=[muAllorig;mean(wholeOrigWm(ind))];
    
    
    
end
SITES{1}.muOrig=mean(muAllorig);

s=2;
muAllorig=[]; muAllharm=[];
for i=1:SITES{s}.ImageNum
    data  = SITES{s}.InputImages{i};
    load([data.Image_Harmonized_dir,'newfile_suffix.mat']);
    load([data.Image_Harmonized_dir,'newmask_suffix.mat']);
    load([data.Image_Harmonized_dir,'processes.mat']);
    SITES{s}.InputImages{i}.DWISuffix = newfile_suffix;
    SITES{s}.InputImages{i}.MaskSuffix = newmask_suffix;
    SITES{s}.InputImages{i}.denoised = process(1);
    SITES{s}.InputImages{i}.bValueMapped = process(2);
    SITES{s}.InputImages{i}.resampled = process(3);
    
    data  = SITES{s}.InputImages{i};
    
    % data.DWISuffix: can be '_rescaled', '_bmapped', or '_denoised'
    if ~isempty(data.DWISuffix)
        data.ImageName = [data.ImageName, data.DWISuffix];
        data.MaskName = [data.MaskName, data.MaskSuffix];
        data.MaskFullPath =  [data. Image_Harmonized_dir, data.MaskName, '.', data.MaskImageType];
        data.ImageFullPath = [data.Image_Harmonized_dir, data.ImageName, '.',  data.ImageType];
    end
    
    
    dataImg = [data.ImageDirectory 'FSL/' SITES{s}.InputImages{i}.ImageName  '_InMNI_FA.nii.gz'];
    harmonizedImg=  [data.ImageDirectory 'FSL/harmonized_' data.ImageName  '_InMNI_FA.nii.gz'];
    
    niiOrig = load_untouch_nii(dataImg);
    niiharmonized = load_untouch_nii(harmonizedImg);
    
    ind=atlas.img == 1;
    wholeOrigWm=niiOrig.img.*single(atlas.img);
    wholeHarmWm=niiharmonized.img.*single(atlas.img);
    
    muAllorig=[muAllorig;mean(wholeOrigWm(ind))];
    muAllharm=[muAllharm;mean(wholeHarmWm(ind))];
    
end
SITES{2}.muOrig=mean(muAllorig);
SITES{2}.muHarm=mean(muAllharm);

end


