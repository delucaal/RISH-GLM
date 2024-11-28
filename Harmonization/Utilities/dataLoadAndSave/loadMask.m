function [mask, sd, moption] = loadMask(maskfullpath, maskimagetype, fid)

   if strcmp(maskimagetype,'nhdr') || strcmp(maskimagetype,'nrrd') 
        mask = nrrdLoadWithMetadata(maskfullpath);
        moption =  mask;
        sd  = mask.spacedirections;
      
         if strcmp(mask.space,'LPS') || strcmpi(mask.space,'left-posterior-superior')
            RAS = [-1 0 0;0 -1 0;0 0 1];
            sd = RAS * sd;
         end
         mask = mask.data;
       
        

    %elseif(strcmp(maskimagetype,'mat'))
    %    load(maskfullpath);
    elseif(strcmp(maskimagetype,'nii') || strcmp(maskimagetype,'img')||...
        strcmp(maskimagetype,'img.gz') || strcmp(maskimagetype,'nii.gz'))       
        mask = load_untouch_nii(maskfullpath);
        moption=mask;
        Sm=diag(mask.hdr.dime.pixdim(2:4));
        sd = Sm(1:3,1:3);
        mask = mask.img;
        
   else
       fprintf(fid,'Warning: no mask');
       mask=ones(size(1), size(2), size(3));
       sd=[];
   end
    
        
end