function flag = verifyData(filepath, filetype, fid)


flag  = tryLoadingDWI(filepath, filetype, fid);


end


function  flag =tryLoadingDWI(imagefullpath, imagetype, fid)
    flag = 1;
    if strcmp(imagetype, 'nhdr') ||  strcmp(imagetype, 'nrrd') 
        try
            dwi = nrrdLoadWithMetadata(imagefullpath); 
        catch
           fprintf(fid, [imagefullpath, '\n']);
           flag = 0;
        end
      
    elseif(strcmp(imagetype, 'nii') || strcmp(imagetype, 'nii.gz')...
            || strcmp(imagetype, 'img') || strcmp(imagetype, 'img.gz')) %single shell only
        
        try
            nii = load_untouch_nii(imagefullpath); 
            
        catch
           fprintf(fid, [imagefullpath, '\n']);
            flag = 0;
        end
      
    end

end