function flag = verifyScaleMap(filepath, filetype, fid)


flag  = tryLoadingScaleMap(filepath, filetype, fid);


end


function  flag =tryLoadingScaleMap(imagefullpath, imagetype, fid)
    flag = 1;
    if strcmp(imagetype, 'nhdr') ||  strcmp(imagetype, 'nrrd') 
        try
            sc = nrrdLoadWithMetadata(imagefullpath); 
            if mean(sc.data(:))==0
                flag=0;
            end
        catch
           fprintf(fid, [imagefullpath, '\n']);
           flag = 0;
        end
      
    elseif(strcmp(imagetype, 'nii') || strcmp(imagetype, 'nii.gz')...
            || strcmp(imagetype, 'img') || strcmp(imagetype, 'img.gz')) %single shell only  
        try
            nii = load_untouch_nii(imagefullpath); 
            if mean(nii.img(:))==0
                flag=0;
            end
        catch
           fprintf(fid, [imagefullpath, '\n']);
            flag = 0;
        end
      
    end

end