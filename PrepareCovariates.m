
load('/autofs/arch11/DATA/PROVIDI_LAB/Alexander/My_data/FTD_RISC/ProcessedData/DEMOGRAPHICS_MATCHED_INC_PRECONSYM_NODUP_WEC_AGEC2_TS.mat');
input_files = {'Source_LEID_1','Source_LEID_3_4','Source_RTM_5_6'};

covariates = [];
caselist = {};
site = [];

for x=1:length(input_files)
   f = fopen([input_files{x} '_dwi.txt'],'rt');
   while(~feof(f))
       l = fgetl(f);
       caselist = cat(1,caselist,{l});
   end
   fclose(f);
end

for x=1:length(caselist)
    [fp,fn,ext] = fileparts(caselist{x});
    for ix=0:6
        fn = strrep(fn,['FU' num2str(ix)],['FU_' num2str(ix)]);
    end
    idx = locateindb(DEMOGRAPHICS_MATCHED, fn);
    sex = strtrim(DEMOGRAPHICS_MATCHED{idx,6});
    if(strcmpi(sex,'man') || strcmpi(sex,'male'))
        sex = 1;
    elseif(strcmpi(sex,'vrouw') || strcmpi(sex,'female'))
        sex = 0;
    else
        error('unexpected');
    end
    covariates = cat(1,covariates,[DEMOGRAPHICS_MATCHED{idx,8} sex]);
    site = cat(1,site,DEMOGRAPHICS_MATCHED{idx,17});
end

% % R5 or R4
% covariates(:,end+1) = site==1; % not needed it is the main constrast
% 8ch
% covariates = [site==3 covariates];
save('GLM_Covariates_Fix.txt','covariates','-ascii')

function idx = locateindb(db,name)
    idx = -1;
    stp = strfind(name,'_dwi_FP');
    std = strfind(name,'sub-');
    name = name(std:stp-1);
    for ix=2:size(db)
       if(strcmp(db{ix,2},name))
           idx = ix;
           return;
       end
    end
end
