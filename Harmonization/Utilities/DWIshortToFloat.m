% dwi short to float conversion
% InputDWIshortnrrd = input filename (short)
% OutputDWIfloatnrrd = output filename (float)
function DWIshortToFloat(InputDWIshortnrrdList, OutputDWIfloatnrrdList)

[~,~,input] = xlsread(InputDWIshortnrrdList);
[~,~,output] = xlsread(OutputDWIfloatnrrdList);

for i=1:size(input,1)
dwi=nrrdLoadWithMetadata(input{i,1});
dwifinal = dwi;
dwifinal.data = single(dwi.data);
nrrdSaveWithMetadata(output{i,1}, dwifinal);
end

end
