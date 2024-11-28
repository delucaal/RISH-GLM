function [MALE, FEMALE, missing] = computeGA(filename )
%filename = '/projects/schiz/suheyla/Data/UPENN_CIDAR_post/After_FA_Phil_test_All.xlsx';

[~,~,raw] = xlsread(filename);



for i=1:size(raw,1)
    
    if  exist(raw{i,1},'file')
        fa = 
        
        FA_all = []; volumeAll=[];
        for b=1:length(majorBundles)
            FA0 =  fa.img.*double(dataBundles(:,:,:,b)>0.0); 
            FA1 =  fa.img.*double(dataBundles(:,:,:,b)>0.1);
            FA2 =  fa.img.*double(dataBundles(:,:,:,b)>0.2);
         
            FA_all = [FA_all ; mean(FA0(FA0>0)) mean(FA1(FA1>0)) mean(FA2(FA2>0))];
           
        end
        if raw{i,4}=='F' %Female
            
            FEMALE.major_bundles = cat(3, FEMALE.major_bundles, FA_all);
         
            FEMALE.id = vertcat(FEMALE.id, raw{i,2});
            FEMALE.age = vertcat(FEMALE.age, raw{i,3});
           
            f  = f+1;  
        else
            
            MALE.major_bundles = cat(3, MALE.major_bundles, FA_all);
          
            MALE.id=vertcat(MALE.id, raw{i,2});
            MALE.age=vertcat(MALE.age, raw{i,3});
            
            m  = m+1;
            
        end
       
    else
        disp([num2str(raw{i,2}) 'is missing']);
        missing.id =  vertcat(missing.id, raw{i,2});
        missing.age =  vertcat(missing.age, raw{i,3});
        missing.gender =  vertcat(missing.gender, raw{i,4});
        mc = mc+1;
    end
    disp(['i: ', num2str(i), ' m: ', num2str(m), ' f: ' num2str(f), ' missing: ', num2str(mc)]);
end
   
   