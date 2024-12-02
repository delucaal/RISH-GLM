%% This function creates RISH feature templates for each site

function createTemplatesUsingRegistration(SITES,templatesPath,option)
SHorder = option.SHOrder;
ANTweights = [1 1];
% % %
% %STEP1: Create a list of RISH feature images for all data

reference_template = option.reference_template;

if ~exist([templatesPath  'Mean_' SITES{1}.name '_FA.nii.gz'], 'file') ||  ~exist([templatesPath  'Mean_' SITES{2}.name '_FA.nii.gz'], 'file')...
        || option.force
    
    fid1=fopen([templatesPath 'CaseList_MultivariateTemplate.txt'],'w');
    for s=1:length(SITES)
        
        for i=1:SITES{s}.ImageNum
            a = [];
            
            SITES{s}.InputImages{i}.ImageName = [SITES{s}.InputImages{i}.ImageName SITES{s}.InputImages{i}.DWISuffix];
            SITES{s}.InputImages{i}.MaskName = [SITES{s}.InputImages{i}.MaskName SITES{s}.InputImages{i}.MaskSuffix];
            
            
            dtiPath=[SITES{s}.InputImages{i}.ImageDirectory, option.DTIdir, SITES{s}.InputImages{i}.ImageName];
            a =horzcat(a, [dtiPath,  '_FA.nii.gz,'] );
            
            
            rishPath=[SITES{s}.InputImages{i}.Image_Harmonized_dir, 'rish/', SITES{s}.InputImages{i}.ImageName];
            
            for l=0:2:(length(ANTweights)-2)*2
                a =horzcat(a, [rishPath, '_L' num2str(l) '.nii.gz,'] );
            end
            fprintf(fid1,'%s\n',a);
            
                
            if(isempty(dir([templatesPath  SITES{s}.InputImages{i}.ImageName '*'])))
                fix=reference_template;%[ templatesPath, 'Mean_', SITES{s}.name '_FA.nii.gz'];
                mov=[SITES{s}.InputImages{i}.ImageDirectory option.DTIdir  SITES{s}.InputImages{i}.ImageName '_FA.nii.gz'];
                [status,out] =system(['antsRegistrationSyNQuick.sh -d 3 -f'...
                    fix  ' -m '  mov  ' -o '  templatesPath  SITES{s}.InputImages{i}.ImageName  '_FA']);

                if status
                    fprintf(fidH, '%s\n', out); fclose(fidH);
                    error(['Error in antsRegistrationSyNQuick %s\n',  data.ImageFullPath]);
                end
            end
        end
    end
    fclose(fid1);
    
    %     % STEP2: Run antsMultivariateTemplateConstruction2.sh
    %     % d: number of dimensions
    %     % i:  Iteration limit (default 4): iterations of the template construction
    %     % k: number of modalities (L0, L2, L4, L6, L8)
    %     % w: Modality weights used in the similarity metric (default = 1)
    %     % t: Type of tranformation model for registration (SyN = Greedy)
    %     % m: Type of similarity metric (CC =  Normalized Cross Correlation)
    %     % o: Output template
    strW = num2str(ANTweights(1));
    for l=2:length(ANTweights)
        strW = horzcat(strW, [ 'x' num2str(ANTweights(l))]);
    end
           
    
    % STEP3: Create a list of *Warp.nii.gz and *GenericAffine.mat for antsApplyTransforms.sh
    currentFolder = pwd;
    cd(templatesPath );
    system('ls *Warp.nii.gz |  grep -v "InverseWarp" > warpPaths.txt');
    system('ls *GenericAffine.mat |  grep -v "template"  > transformPaths.txt');
    cd(currentFolder);
    
    fidW = fopen([templatesPath 'warpPaths.txt'],'r');
    fidT = fopen([templatesPath 'transformPaths.txt']);
    wline = fgets(fidW); tline = fgets(fidT);
    c=1;
    while ischar(wline) && ischar(tline)
        warps{c}=wline(1:end-1); transforms{c}=tline(1:end-1);
        tline = fgets(fidT); wline = fgets(fidW);
        c=c+1;
    end
    fclose(fidW); fclose(fidT);
    if(length(warps) ~= length(transforms))
        fprintf(['Something went wrong while creating templates.' ...
            'Please delete your template folder and re-run the method.']);
        return;
    end
    
    system(['fslmerge -t ' templatesPath '/all_warped_FA ' templatesPath '/*FAWarped.nii.gz']);
    system(['fslmaths ' templatesPath '/all_warped_FA -Tmean ' templatesPath '/template0']);
    % % %
    % % %
    % % % %
    % % % % % %    STEP 4: Warp all the bands: L0 L2 L4 L6 L8 ...
    diffusionMeasures = { 'MD', 'FA'};
    for s=1:length(SITES)
        for i=1:SITES{s}.ImageNum
            
            rishPath=[SITES{s}.InputImages{i}.Image_Harmonized_dir, 'rish/', SITES{s}.InputImages{i}.ImageName];
            dtiPath=[SITES{s}.InputImages{i}.ImageDirectory,option.DTIdir, SITES{s}.InputImages{i}.ImageName];
            
            if ~isempty(SITES{s}.InputImages{i}.MaskSuffix)
                movMaskPath = [SITES{s}.InputImages{i}.Image_Harmonized_dir,SITES{s}.InputImages{i}.MaskName, '.', SITES{s}.InputImages{i}.MaskImageType];
            else
                movMaskPath = [SITES{s}.InputImages{i}.ImageDirectory,SITES{s}.InputImages{i}.MaskName,  '.', SITES{s}.InputImages{i}.MaskImageType];
            end
            warpedMaskPath = [SITES{s}.InputImages{i}.Image_Harmonized_dir,SITES{s}.InputImages{i}.MaskName, 'Warped.nii.gz'];
            
            imgName = [SITES{s}.InputImages{i}.ImageName, '_FA'];
            ispresent = cellfun(@(st) ~isempty(strfind(st,imgName)), warps); indW=find(ispresent);
            ispresent = cellfun(@(st) ~isempty(strfind(st,imgName)), transforms); indT=find(ispresent);
            
            if(~isempty(indW) && ~isempty(indT))
                
                for mi=1:numel(diffusionMeasures)
                    system(['antsApplyTransforms -d 3 -i ' dtiPath, '_', diffusionMeasures{mi}, '.nii.gz'...
                        ' -o '  dtiPath, '_Warped', diffusionMeasures{mi}, '.nii.gz' ...
                        ' -r ' templatesPath,   'template0.nii.gz'...
                        ' -t ' templatesPath,  warps{indW} ...
                        ' -t ' templatesPath,  transforms{indT}]);
                    
                    
                end
                system(['antsApplyTransforms -d 3 -i ' movMaskPath...
                    ' -o '  warpedMaskPath...
                    ' -r ' templatesPath,   'template0.nii.gz'...
                    ' -t ' templatesPath,  warps{indW} ...
                    ' -t ' templatesPath,  transforms{indT}]);
                for l=0:2:SHorder
                    
                    system(['antsApplyTransforms -d 3 -i ' rishPath, '_L' num2str(l) '.nii.gz'...
                        ' -o '  rishPath, '_WarpedL' num2str(l) '.nii.gz'...
                        ' -r ' templatesPath,   'template0.nii.gz'...
                        ' -t ' templatesPath,  warps{indW} ...
                        ' -t ' templatesPath,  transforms{indT}]);
                end
            end
        end
    end
    
    
    %
    % %   STEP 5: Load the warped images and compute the templates at different sites/orders
    %
    %
    % %% Mean diffusion measures
    for s=1:length(SITES)
        
        
        for mi=1:numel(diffusionMeasures)
            temp{1}=[];
            mask{1}=[];
            
            for i=1:SITES{s}.ImageNum
                
                dtiPath=[SITES{s}.InputImages{i}.ImageDirectory,option.DTIdir, SITES{s}.InputImages{i}.ImageName];
                disp(dtiPath);
                nii=load_untouch_nii([dtiPath, '_Warped', diffusionMeasures{mi}, '.nii.gz']);
                temp{1}= cat(4,temp{1}, (nii.img));
                if mi==1
                    niiMask=load_untouch_nii([SITES{s}.InputImages{i}.Image_Harmonized_dir,SITES{s}.InputImages{i}.MaskName, 'Warped.nii.gz']);
                    mask{1}= cat(4,mask{1}, (niiMask.img));
                end
            end
            
            nii.img = mean(temp{1},4);
            save_untouch_nii(nii,[templatesPath  'Mean_' SITES{s}.name '_', diffusionMeasures{mi},  '.nii.gz']);
            
            nii.img = var(temp{1},1,4);
            save_untouch_nii(nii,[templatesPath  'Std_' SITES{s}.name '_', diffusionMeasures{mi},  '.nii.gz']);
            
            if mi==1
                nii.img = single(mean(mask{1},4)>0.5);
                se = strel3d(3);
                nii.img = imopen(nii.img, se);
                save_untouch_nii(nii,[templatesPath   SITES{s}.name '_Mask.nii.gz']);
            end
            
            
            meanDiff = temp{1};
            save([templatesPath  SITES{s}.name '_', diffusionMeasures{mi}, '.mat'], 'meanDiff');
        end
        
    end
    
    maskRef = load_untouch_nii([templatesPath   SITES{1}.name '_Mask.nii.gz']);
    maskOther = load_untouch_nii([templatesPath   SITES{2}.name '_Mask.nii.gz']);
    
    %mask
    maskRef.img = maskRef.img.*maskOther.img;
    save_untouch_nii(maskRef,[templatesPath   'templateMask.nii.gz']);
    for s=1:length(SITES)
        for mi=1:numel(diffusionMeasures)
            nii = load_untouch_nii([templatesPath  'Mean_' SITES{s}.name '_', diffusionMeasures{mi},  '.nii.gz']);
            nii.img = nii.img.*maskRef.img;
            save_untouch_nii(nii,[templatesPath  'Mean_' SITES{s}.name '_', diffusionMeasures{mi},  '.nii.gz']);
            
            nii = load_untouch_nii([templatesPath  'Std_' SITES{s}.name '_', diffusionMeasures{mi},  '.nii.gz']);
            nii.img = nii.img.*maskRef.img;
            save_untouch_nii(nii,[templatesPath  'Std_' SITES{s}.name '_', diffusionMeasures{mi},  '.nii.gz']);
        end
    end
    
    
    
    % %% Diffusion Measures (GFA, GFA, MD) differences:
    %
    % For traveling heads,  subtract other site from reference site for each
    % subject. Then take the mean.
    for mi=1:numel(diffusionMeasures)
        if option.travelingHeads
            temp{1}=[]; temp{2}=[]; temp{3}=[];
            
            for i=1:SITES{1}.ImageNum
                
                dtiPathRef=[SITES{1}.InputImages{i}.ImageDirectory,option.DTIdir, SITES{1}.InputImages{i}.ImageName];
                niiRef=load_untouch_nii([dtiPathRef, '_Warped', diffusionMeasures{mi}, '.nii.gz']);
                
                dtiPathOther=[SITES{2}.InputImages{i}.ImageDirectory,option.DTIdir, SITES{2}.InputImages{i}.ImageName];
                niiOther=load_untouch_nii([dtiPathOther, '_Warped', diffusionMeasures{mi}, '.nii.gz']);
                
                
                temp{1}= cat(4,temp{1}, ( niiRef.img-niiOther.img).*single(maskRef.img));
                
                perct = 100*(niiRef.img-niiOther.img)./(niiRef.img+eps);
                perct(isnan(perct))=0; perct(perct>100)=100; perct(perct<-100)=-100;
                temp{2}= cat(4,temp{2}, perct.*single(maskRef.img));
                
                perct = 100*(smooth3(niiRef.img)-smooth3(niiOther.img))./(smooth3(niiRef.img)+eps);
                perct(isnan(perct))=0; perct(perct>100)=100; perct(perct<-100)=-100;
                temp{3}= cat(4,temp{3}, perct.*single(maskRef.img));
                
            end
            niiOther.img = mean(temp{1},4);
            save_untouch_nii(niiOther,[templatesPath  'Delta_', diffusionMeasures{mi}, '.nii.gz']);
            
            mu = mean(temp{2},4);
            mu(mu>100)=100; mu(isnan(mu))=0;
            niiOther.img = mu;
            save_untouch_nii(niiOther,[templatesPath  'PercentageDiff_', diffusionMeasures{mi}, '.nii.gz']);
            clear mu;
            
            mu = mean(temp{3},4);
            mu(mu>100)=100; mu(isnan(mu))=0;
            niiOther.img =  mu;
            save_untouch_nii(niiOther,[templatesPath  'PercentageDiff_', diffusionMeasures{mi}, 'smooth.nii.gz']);
            
        else
            % For others, subtract the mean of other site from the mean of
            % reference site.
            niiRef = load_untouch_nii([templatesPath  'Mean_' SITES{1}.name '_', diffusionMeasures{mi}, '.nii.gz']);
            niiOther = load_untouch_nii([templatesPath  'Mean_' SITES{2}.name '_', diffusionMeasures{mi}, '.nii.gz']);
            
            niiRef.img = niiRef.img.*single(maskRef.img);
            niiOther.img = niiOther.img.*single(maskRef.img);
            nii = niiOther;
            nii.img = (niiRef.img-niiOther.img);
            save_untouch_nii(nii,[templatesPath  'Delta_', diffusionMeasures{mi}, '.nii.gz']);
            
            perct = 100*(niiRef.img-niiOther.img)./(niiRef.img+eps);
            perct(isnan(perct))=0; perct(perct>100)=100; perct(perct<-100)=-100;
            nii.img = perct;
            save_untouch_nii(nii,[templatesPath  'PercentageDiff_', diffusionMeasures{mi}, '.nii.gz']);
            
            perct = 100*(smooth3(niiRef.img)-smooth3(niiOther.img))./(smooth3(niiRef.img)+eps);
            perct(isnan(perct))=0; perct(perct>100)=100; perct(perct<-100)=-100;
            nii.img =  perct;
            save_untouch_nii(nii,[templatesPath  'PercentageDiff_', diffusionMeasures{mi}, 'smooth.nii.gz']);
            
        end
    end
    
    % %% Mean Rish Templates
    for s=1:length(SITES)
        for l=0:2:SHorder
            temp{1}=[];
            clear nii;
            for i=1:SITES{s}.ImageNum
                rishPath=[SITES{s}.InputImages{i}.Image_Harmonized_dir, 'rish/', SITES{s}.InputImages{i}.ImageName];
                nii=load_untouch_nii([rishPath, '_WarpedL' num2str(l) '.nii.gz']);
                temp{1}= cat(4,temp{1}, (nii.img));
            end
            nii.img = mean(temp{1},4);
            save_untouch_nii(nii,[templatesPath  'Mean_' SITES{s}.name '_L'  num2str(l)  '.nii.gz']);
            
            nii.img = var(temp{1},1, 4);
            save_untouch_nii(nii,[templatesPath  'Std_' SITES{s}.name '_L'  num2str(l)  '.nii.gz']);
            
        end
    end
    
    % %% Difference Rish
    %
    % % For traveling heads,  subtract other site from reference site for each
    % % subject. Then take the mean.
    
    
    if option.travelingHeads
        
        for l=0:2:SHorder
            delta=[]; scale=[];perct=[];  percts=[];
            for i=1:SITES{1}.ImageNum
                refPath=[SITES{1}.InputImages{i}.Image_Harmonized_dir, 'rish/', SITES{1}.InputImages{i}.ImageName];
                refImg=load_untouch_nii([refPath, '_WarpedL' num2str(l) '.nii.gz']);
                
                movPath=[SITES{2}.InputImages{i}.Image_Harmonized_dir, 'rish/', SITES{2}.InputImages{i}.ImageName];
                movImg=load_untouch_nii([movPath, '_WarpedL' num2str(l) '.nii.gz']);
                
                delta= cat(4,delta,((refImg.img)-(movImg.img)).*single(maskRef.img));
                scale= cat(4,scale,((refImg.img)./(movImg.img+eps)).*single(maskRef.img));
                
                diff = 100*(refImg.img-movImg.img)./(refImg.img+eps);
                diff(isnan(diff))=0; diff(diff>100)=100; diff(diff<-100)=-100;
                perct= cat(4,perct,diff.*single(maskRef.img));
                
                diff = 100*(smooth3(refImg.img)-smooth3(movImg.img))./(smooth3(refImg.img)+eps);
                diff(isnan(diff))=0; diff(diff>100)=100; diff(diff<-100)=-100;
                percts= cat(4,percts,diff.*single(maskRef.img));
                
            end
            
            refImg.img = sqrt(mean(scale,4));
            save_untouch_nii(refImg,[templatesPath  'Scale_L'  num2str(l) '.nii.gz']);
            
            movImg.img = mean(delta,4);
            save_untouch_nii(movImg,[templatesPath  'Delta_L'  num2str(l) '.nii.gz']);
            
            refImg.img = mean(perct, 4);
            save_untouch_nii(refImg,[templatesPath  'PercentageDiff_L'  num2str(l) '.nii.gz']);
            
            refImg.img =  mean(percts, 4);
            save_untouch_nii(refImg,[templatesPath  'PercentageDiff_L'  num2str(l) 'smooth.nii.gz']);
        end
        
    else
        % For others, subtract the mean of other site from the mean of
        % reference site.
        for l=0:2:SHorder
            
            
            refImg=load_untouch_nii([templatesPath  'Mean_' SITES{1}.name '_L'  num2str(l)  '.nii.gz']);
            movImg=load_untouch_nii([templatesPath  'Mean_' SITES{2}.name '_L'  num2str(l)  '.nii.gz']);
            
            
            refImg.img = refImg.img.*single(maskRef.img);
            movImg.img = movImg.img.*single(maskRef.img);
            
            delta1 = refImg.img - movImg.img;
            scale1 = sqrt(refImg.img./(movImg.img+eps));
            refTmp = refImg.img; movTmp = movImg.img;
            
            refImg.img = scale1.*single(maskRef.img);
            movImg.img = delta1.*single(maskRef.img);
            save_untouch_nii(refImg,[templatesPath  'Scale_L'  num2str(l) '.nii.gz']);
            save_untouch_nii(movImg,[templatesPath  'Delta_L'  num2str(l) '.nii.gz']);
            
            
            mu = 100*(refTmp - movTmp)./(refTmp+eps);
            mu(mu>100)=100; mu(isnan(mu))=0; mu(mu<-100)=-100;
            refImg.img = mu.*single(maskRef.img);
            save_untouch_nii(refImg,[templatesPath  'PercentageDiff_L'  num2str(l) '.nii.gz']);
            
            mu = 100*(smooth3(refTmp) - smooth3(movTmp))./(smooth3(refTmp)+eps);
            mu(mu>100)=100; mu(isnan(mu))=0;mu(mu<-100)=-100;
            refImg.img = mu.*single(maskRef.img);
            save_untouch_nii(refImg,[templatesPath  'PercentageDiff_L'  num2str(l) 'smooth.nii.gz']);
            clear mu;
        end
        
        
        
    end
    
end


end
