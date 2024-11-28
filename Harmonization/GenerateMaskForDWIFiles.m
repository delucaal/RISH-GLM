function GenerateMaskForDWIFiles(folder)
    call = 'bet';
    content = dir(fullfile(folder,'*trafo_rek.nii.gz'));
    
    for ij=1:length(content)
       fname = fullfile(content(ij).folder,content(ij).name(1:end-7));
       if(exist([fname '_brain_mask.nii.gz'],'file') < 1)
          cmd = [call ' ' fname ' ' fname '_brain  -m -f 0.3 -R '];
          final_cmd = bashify_cmd(cmd);
          system(final_cmd);
       end
    end
end

function a = PathsLinux2Windows(a)
    windows_drive = 'y:';
    windows_drive_cap = upper(windows_drive);
    linux_drive = '/mnt/y/';
    a = strrep(a, windows_drive, linux_drive);
    a = strrep(a, windows_drive_cap, linux_drive);
    a = strrep(a,'\','/');
end

function a = PathsWindows2Linux(a)
    windows_drive = 'z:';
    linux_drive = '/mnt/z/';
    a = strrep(a, linux_drive, windows_drive);
    a = strrep(a,'/','\');
end

function cmd = bashify_cmd(cmd)
    if(isunix)
        return
    end
    cmd = ['bash -c "' cmd '"'];
    cmd = PathsLinux2Windows(cmd);
end
