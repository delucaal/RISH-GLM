%% This function reads _dwi.txt and _mask.txt files (it should be prepared seperately for each site by user)
%% Each file should contain full path information of dMRIs and brain masks.
function SITES =loadData(SITES, option)


fid1 = fopen(['./log/log_', option.name,'loadData1'],'w');
if (isempty(SITES))
    fprintf(fid1,'No site is entered');
    error(['./log/log_', option.name,'loadData1: No site is entered']);
end

fidV = fopen(['./log/log_', option.name,'verifyData1'],'w');

%disp('Loading and verifying data...');
for s=1:length(SITES)
    %disp(SITES{s}.fullpath);
    %[~,~,dwipaths] = xlsread(SITES{s}.fullpath);
    dwipaths = readTXTFiles(SITES{s}.fullpath);
    if(isempty(dwipaths))
        fprintf(fid1,'Empty file');
        return;
    else
        t = size(dwipaths);
        
        if(t(2)<2)
            fprintf(fid1,'Input information is missing');
            return;
        elseif(t(2)==2)
            SITES{s}.ImageNum = t(1);
            for i=1:t(1)
                
                if exist(dwipaths{i,1}, 'file') &&  exist(dwipaths{i,2}, 'file')
                    idxSlash = strfind(dwipaths{i,1},'/'); imgpath = dwipaths{i,1}; %find the image file name, e.g, 'Id.nii.gz'
                    imageFile = imgpath(idxSlash(end)+1:end);
                    idxDot=strfind(imageFile,'.'); %find the image name, e.g, 'Id'
                    SITES{s}.InputImages{i}.ImageName=imageFile(1:idxDot(1)-1);
                    SITES{s}.InputImages{i}.ImageType = imageFile(idxDot(1)+1:end);
                    SITES{s}.InputImages{i}.ImageDirectory=imgpath(1:idxSlash(end));
                    SITES{s}.InputImages{i}.ImageFullPath=imgpath;
                    
                    idxSlash = strfind(dwipaths{i,2},'/'); maskpath=dwipaths{i,2};
                    masktype = maskpath(idxSlash(end)+1:end);
                    idxDot=strfind(masktype,'.');
                    SITES{s}.InputImages{i}.MaskName=masktype(1:idxDot(1)-1);
                    SITES{s}.InputImages{i}.MaskFullPath=maskpath;
                    SITES{s}.InputImages{i}.MaskDirectory=maskpath(1:idxSlash(end));
                    SITES{s}.InputImages{i}.MaskImageType=masktype(idxDot(1)+1:end);
                   
                    SITES{s}.InputImages{i}.resampled = 0;
                    SITES{s}.InputImages{i}.bValueMapped = 0;
                    SITES{s}.InputImages{i}.denoised = 0;
                    SITES{s}.InputImages{i}.DWISuffix = '';
                    SITES{s}.InputImages{i}.MaskSuffix = '';
                    SITES{s}.InputImages{i}.harmonized = 0;
                    
                   %verifyData( SITES{s}.InputImages{i}.ImageFullPath, SITES{s}.InputImages{i}.ImageType, fidV);
                   
                    if ~isempty(option.harmonizeDir)
                        harmpath = option.harmonizeDir;
                        
                    else
                        
                        cnt=0; harmpath=['harmonized' num2str(cnt) '_'];
                        while  exist([SITES{s}.InputImages{i}.ImageDirectory harmpath '/' ],'file') && ~option.force
                            cnt=cnt+1;
                            harmpath=['harmonized' num2str(cnt) '_'];
                        end
                    end
                    if ~exist([SITES{s}.InputImages{i}.ImageDirectory harmpath '/'], 'dir')
                        mkdir([SITES{s}.InputImages{i}.ImageDirectory harmpath '/']);
                    end
                    SITES{s}.InputImages{i}.Image_Harmonized_dir=[SITES{s}.InputImages{i}.ImageDirectory harmpath '/'];
                   % i
                    fprintf(fid1,'%s\n', SITES{s}.InputImages{i}.ImageFullPath);
                else
                     fprintf(fid1,'%s or %s does not exist\n', dwipaths{i,1}, dwipaths{i,2});
                     error(['error in loading data, check ', './log/log_', option.name,'loadData1' ' for details']);
                end
            end
        else
            fprintf('Wrong additional column');
        end
    end
    %s
    
end

fclose(fid1);
fclose(fidV);
dirVer = dir(['./log/log_', option.name,'verifyData1']);
if dirVer.bytes > 0
     error(['check ./log/log_', option.name,'verifyData1']);
end
end


function paths = readTXTFiles(filefullpath)
fid = fopen([filefullpath, '_dwi.txt']);
tline =  fgetl(fid);
i = 1;
while ~feof(fid) && ~isempty(tline) && ~strcmpi(tline, 'nan')
    paths{i,1} = tline;
    i =  i+1;
    tline =  fgetl(fid);
end
if  ~isempty(tline) && ~sum(isspace(tline)) &&  ~strcmpi(tline, 'nan')
    paths{i,1} = tline;
end
fclose(fid);

fid = fopen([filefullpath, '_mask.txt']);
tline =  fgetl(fid);
i = 1;
while ~feof(fid)  && ~isempty(tline) && ~strcmpi(tline, 'nan')
    paths{i,2} = tline;
    i =  i+1;
    tline =  fgetl(fid);
end
if  ~isempty(tline) &&   ~sum(isspace(tline)) &&  ~strcmpi(tline, 'nan')
    paths{i,2} = tline;
end
fclose(fid);
end















