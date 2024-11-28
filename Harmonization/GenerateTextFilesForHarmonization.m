function GenerateTextFilesForHarmonization(directory,out_file_prefix,max_elements,filter)
    if(exist('filter','var') < 1)
        filter = '';
    end

    files_dwi = dir(fullfile(directory,['*' filter '*_trafo.nii*']));
    files_mask = dir(fullfile(directory,['*' filter '*_mask.nii*']));

    if(isempty(files_dwi) || length(files_dwi) ~= length(files_mask))
        error('Something is wrong in the indicated directory');
    end
    
    if(nargin < 3)
        max_elements = length(files_dwi);
    end
        
    if(isfinite(max_elements))
        step = ceil(length(files_dwi)/max_elements);
    else
        step = 1;
    end
    
    files_dwi = files_dwi(1:step:end);
    files_mask = files_mask(1:step:end);
    
    if(length(files_dwi) > max_elements)
        files_dwi = files_dwi(1:max_elements);
        files_mask = files_mask(1:max_elements);
    end
    
    file_dwi = fopen([out_file_prefix '_dwi.txt'],'wt');
    file_mask = fopen([out_file_prefix '_mask.txt'],'wt');
    
    for ij=1:min(length(files_dwi),max_elements)
        files_dwi(ij).folder = directory; %bruno
        files_mask(ij).folder = directory; %bruno changed - dir is not reading the folder of the file and 'newline' does not work

       fprintf(file_dwi,'%s%s%s',fullfile(files_dwi(ij).folder,files_dwi(ij).name));
       fprintf(file_mask,'%s%s%s',fullfile(files_mask(ij).folder,files_mask(ij).name));
    end
    
    fclose(file_dwi);
    fclose(file_mask);
    
end