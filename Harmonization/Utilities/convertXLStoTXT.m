function convertXLStoTXT(pathxls)
[~,~,raw] = xlsread(pathxls);
[folder, name, ~] = splitFilePath(pathxls);

fid = fopen([folder, name, '_dwi.txt'], 'w');

for i=1:size(raw,1)
    fprintf(fid, '%s\n', raw{i,1});
end
fclose(fid);

fid = fopen([folder, name, '_mask.txt'], 'w');
for i=1:size(raw,1)
    fprintf(fid, '%s\n', raw{i,2});
end
fclose(fid);

end


% splitting  filepath to substrings: folder, dataname, type
function [folder, name, type] =  splitFilePath(filename)

indt = strfind(filename, '.');
type = filename(indt(1)+1:end);

indf = strfind(filename, '/');
folder = filename(1:indf(end));

name = filename(indf(end)+1:indt(1)-1);


end