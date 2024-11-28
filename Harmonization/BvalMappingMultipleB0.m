input_folder = 'InputFolder';
output_folder = 'OutputFolder';

target_bval = 1000; % the b-value in s/mm^2

target_files = dir(fullfile(input_folder,'*.bval'));
for fid=1:length(target_files)
     bval = load(fullfile(input_folder,target_files(fid).name));
     b0 = bval == 0;
     try
         data = load_untouch_nii(fullfile(input_folder,[target_files(fid).name(1:end-4) 'nii']));
     catch
         data = load_untouch_nii(fullfile(input_folder,[target_files(fid).name(1:end-4) 'nii.gz']));
     end
     data.img = single(data.img);

     dataO = data;
     for bid=1:size(data.img,4)
         if(b0(bid) == 0)
             data.img(:,:,:,bid) = data.img(:,:,:,bid).^(target_bval/bval(bid));
         end
     end

     bval(~b0) = target_bval;
     save_untouch_nii(data,fullfile(output_folder,[target_files(fid).name(1:end-4) 'nii']));
     save(fullfile(output_folder,target_files(fid).name),'bval','-ascii')
     copyfile(fullfile(input_folder,[target_files(fid).name(1:end-4) 'bvec']),fullfile(output_folder,[target_files(fid).name(1:end-4) 'bvec']))
     mask = dir(fullfile(input_folder,[target_files(fid).name(1:end-5) '_mask*']));
     if(~isempty(mask))
         for fid2=1:length(mask)
             copyfile(fullfile(input_folder,mask(fid2).name),fullfile(output_folder,mask(fid2).name))
         end
     end
end
