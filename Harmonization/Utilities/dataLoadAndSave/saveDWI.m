% save DWI data
% DWI can be in three formats: {'nhdr'}, {'mat'} or
% {'nii','img','nii.gz','img.gz'}

%% Inputs:
% imagefullpath : full path the original DWI data file
% harmonized_imagefullpath: full path the new_dwi DWI data file
% imagetype: image type of DWI data file
% new_dwi: struct type of new_dwi DWI data

function saveDWI(imagefullpath, new_dir, new_imgname, imagetype, new_dwi)

    harmonized_imagefullpath = [new_dir new_imgname '.' imagetype];
    
    if(strcmp(imagetype, 'nhdr'))
        dwiHarmonized=struct('gradients', new_dwi.u,'bvalue', max(new_dwi.b),...
            'data',new_dwi.S,'spacedirections',new_dwi.spacedirections,'spaceorigin',new_dwi.spaceorigin);
        
        %mat2DWInhdr(new_imgname, new_dir, dwiHarmonized, 'float'); clear dwi;

    elseif(strcmp(imagetype, 'mat'))
        S_harmonized = new_dwi.S; u = new_dwi.u; b = new_dwi.b;
        s0 = new_dwi.s0; Sm = new_dwi.Sm; voxel = new_dwi.voxel;
        save(harmonized_imagefullpath, 'S_harmonized', 'u', 'b', 's0', 'Sm', 'voxel' );

    elseif(strcmp(imagetype, 'nii') || strcmp(imagetype, 'nii.gz')...
            || strcmp(imagetype, 'img') || strcmp(imagetype, 'img.gz'))
        nii = load_untouch_nii(imagefullpath);
        nii.img = new_dwi.S;
        nii.hdr.dime.dim(5) = size(nii.img,4);
        save_untouch_nii(nii, harmonized_imagefullpath);

        upath = [harmonized_imagefullpath(1:end-length(imagetype)) 'bvec'];
        
        fid = fopen(upath,'w');
        for i=1:length(new_dwi.u)
        fprintf(fid, '%f %f %f\n', new_dwi.u(i,1), new_dwi.u(i,2), new_dwi.u(i,3)); 
        end
        fclose(fid);

        bpath = [harmonized_imagefullpath(1:end-length(imagetype)) 'bval'];
        fid = fopen(bpath,'w'); fprintf(fid, '%f\n', new_dwi.b); fclose(fid);
    else
        fprintf('saveDWI: Wrong file type!');
    end

end