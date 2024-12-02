%% This function creates RISH feature templates for each site

function createTemplateBasedOnExisting_aluca_GLM(SITES,templatesPath,existing_template,option)
SHorder = option.SHOrder;
ANTweights = [1 1];
% % %
% %STEP1: Create a list of RISH feature images for all data
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
    
%         
%     if option.parallel
%         system(['antsMultivariateTemplateConstruction2.sh  -d 3   -i 4  -g 0.2 -k ' num2str(length(ANTweights)) ' -w ', strW  ...
%             ' -t BSplineSyN[0.1,26,0]   -m CC  -r 1 -c 2 -j ', num2str(option.numcores), ' -l 1 -f 8x4x2x1 -s 3x2x1x0 -o ', templatesPath,...
%             ' ', templatesPath, 'CaseList_MultivariateTemplate.txt']);
%     else
%         system(['antsMultivariateTemplateConstruction2.sh  -d 3   -i 4  -g 0.2 -k ' num2str(length(ANTweights)) ' -w ', strW  ...
%             ' -t BSplineSyN[0.1,26,0]   -m CC  -r 1 -c 0 -l 1 -f 8x4x2x1 -s 3x2x1x0 -o ', templatesPath,...
%             ' ', templatesPath, 'CaseList_MultivariateTemplate.txt']);
%     end
    
    system(['ln -s ' existing_template '/*Warp.nii.gz ' templatesPath]);
    system(['ln -s ' existing_template '/*Affine.mat ' templatesPath]);   
    system(['ln -s ' existing_template '/template*.nii.gz ' templatesPath]);
    system(['ln -s ' existing_template '/*Mask*.nii.gz ' templatesPath]);
    system(['ln -s ' existing_template '/*MD*.nii.gz ' templatesPath]);
    system(['ln -s ' existing_template '/*FA*.nii.gz ' templatesPath]);
    
    % STEP3: Create a list of *Warp.nii.gz and *GenericAffine.mat for antsApplyTransforms.sh
    currentFolder = pwd;
    cd(templatesPath );
    system('ls *1Warp.nii.gz > warpPaths.txt');
    system('ls *0GenericAffine.mat > transformPaths.txt');
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
    
    maskRef = load_untouch_nii([templatesPath   SITES{1}.name '_Mask.nii.gz']);
    maskOther = load_untouch_nii([templatesPath   SITES{2}.name '_Mask.nii.gz']);
    
    %mask
    maskRef.img = maskRef.img.*maskOther.img;
    save_untouch_nii(maskRef,[templatesPath   'templateMask.nii.gz']);
    
    % %% Mean diffusion measures
    for mi=1:numel(diffusionMeasures)
        for s=1:length(SITES)
        
        
            temp{s}=[];
            mask{s}=[];
            
            for i=1:SITES{s}.ImageNum
                dtiPath=[SITES{s}.InputImages{i}.ImageDirectory,option.DTIdir, SITES{s}.InputImages{i}.ImageName];
                disp(dtiPath);
                try
                    nii=load_untouch_nii([dtiPath, '_Warped', diffusionMeasures{mi}, '.nii.gz']);
                catch
                    disp('Debug');
                end
                temp{s}= cat(4,temp{s}, (nii.img));
                if mi==1
                    niiMask=load_untouch_nii([SITES{s}.InputImages{i}.Image_Harmonized_dir,SITES{s}.InputImages{i}.MaskName, 'Warped.nii.gz']);
                    mask{s}= cat(4,mask{s}, (niiMask.img));
                end
            end
            
            nii.img = mean(temp{s},4);
            save_untouch_nii(nii,[templatesPath  'Mean_' SITES{s}.name '_', diffusionMeasures{mi},  '.nii.gz']);
            
            nii.img = var(temp{s},1,4);
            save_untouch_nii(nii,[templatesPath  'Var_' SITES{s}.name '_', diffusionMeasures{mi},  '.nii.gz']);
            
            if mi==1
                nii.img = single(mean(mask{s},4)>0.5);
                se = strel3d(3);
                nii.img = imopen(nii.img, se);
                save_untouch_nii(nii,[templatesPath   SITES{s}.name '_Mask.nii.gz']);
            end
            
            
            meanDiff = temp{s};
            save([templatesPath  SITES{s}.name '_', diffusionMeasures{mi}, '.mat'], 'meanDiff');
        end

        covariates = option.covariates;%load('DEMOGRAPHICS_STATS.txt');
        NGroups = length(temp);
        [sx,sy,sz,st] = size(temp{1});
        GroupsLength = cellfun(@(x)size(x,4),temp);
        NPoints = sum(GroupsLength);
        X = zeros(NPoints,NGroups);
        idx = 1;
        for x=1:NGroups
            indices = idx:sum(GroupsLength(1:x));
            idx = idx + GroupsLength(x);
            X(indices,x) = 1;
        end

        V = temp;
        for x=1:length(V)
            V(x) = {reshape(V{x},sx*sy*sz,size(V{x},4))};
        end

        X = cat(2,X,covariates);
        regression = zeros(sx*sy*sz,size(X,2));
        regression_pval = zeros(sx*sy*sz,size(X,2));
        MASK = reshape(maskRef.img,sx*sy*sz,1);
        
        for voxid=1:size(regression,1)
            if(MASK(voxid) > 0)
%                model = fitglm(X,[V1(voxid,:)';V2(voxid,:)'],'intercept',false);
               outliers = [];
               Vl = [];
               for g=1:NGroups
                   if(option.robust_std > 0)
                        [~,o1] = RobustAverage1D(V{g}(voxid,:)',option.robust_std);
                   else
                        o1 = zeros(size(V{g},2),1);
                   end
                   outliers = cat(1,outliers,o1);
                   Vl = cat(2,Vl,V{g}(voxid,:));
               end
%                [b,~,stats] = glmfit(X(~outliers,:),Vl(~outliers),'normal','constant','off');
%                [b,~,stats] = glmfit(X,Vl,'normal','constant','off');
               b = X(~outliers,:)\Vl(~outliers)';
               regression(voxid,:) = b;
%                regression_pval(voxid,:) = stats.p;
            end
        end
        
        nii.img = reshape(regression,sx,sy,sz,size(regression,2));

        if(option.smooth_std > 0)
            for dim4d=1:size(nii.img,4)
                V = nii.img(:,:,:,dim4);
                V = smooth3(V,'gaussian',option.smooth_std);
                nii.img(:,:,:,dim4) = V;
            end
        end
        
        nii.hdr.dime.dim(1) = 4;
        nii.hdr.dime.dim(5) = size(nii.img,4);
        save_untouch_nii(nii,[templatesPath   diffusionMeasures{mi} '_RISHGLM_beta.nii.gz']);
        
        nii.img = reshape(regression_pval,sx,sy,sz,size(regression_pval,2));
        save_untouch_nii(nii,[templatesPath   diffusionMeasures{mi} '_RISHGLM_pval.nii.gz']);
    end
    
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
    for l=0:2:SHorder
        for s=1:length(SITES)
            temp{s}=[];
            clear nii;
            for i=1:SITES{s}.ImageNum
                rishPath=[SITES{s}.InputImages{i}.Image_Harmonized_dir, 'rish/', SITES{s}.InputImages{i}.ImageName];
                nii=load_untouch_nii([rishPath, '_WarpedL' num2str(l) '.nii.gz']);
                temp{s}= cat(4,temp{s}, (nii.img));
            end
            if(option.robust_std > 0)
                nii.img = RobustAverage(temp{s},option.robust_std);
            else
                nii.img = mean(temp{1},4);
            end
            save_untouch_nii(nii,[templatesPath  'Mean_' SITES{s}.name '_L'  num2str(l)  '.nii.gz']);

            nii.img = var(temp{s},1, 4);
            save_untouch_nii(nii,[templatesPath  'Var_' SITES{s}.name '_L'  num2str(l)  '.nii.gz']);
        end
        
        covariates = option.covariates;%load('DEMOGRAPHICS_STATS.txt');

        NGroups = length(temp);
        [sx,sy,sz,st] = size(temp{1});
        GroupsLength = cellfun(@(x)size(x,4),temp);
        NPoints = sum(GroupsLength);
        X = zeros(NPoints,NGroups);
        idx = 1;
        for x=1:NGroups
            indices = idx:sum(GroupsLength(1:x));
            idx = idx + GroupsLength(x);
            X(indices,x) = 1;
        end

        V = temp;
        for x=1:length(V)
            V(x) = {reshape(V{x},sx*sy*sz,size(V{x},4))};
        end

        X = cat(2,X,covariates);
        regression = zeros(sx*sy*sz,size(X,2));
        regression_pval = zeros(sx*sy*sz,size(X,2));
        MASK = reshape(maskRef.img,sx*sy*sz,1);
        
        for voxid=1:size(regression,1)
            if(MASK(voxid) > 0)
%                model = fitglm(X,[V1(voxid,:)';V2(voxid,:)'],'intercept',false);
               outliers = [];
               Vl = [];
               for g=1:NGroups
                   if(option.robust_std > 0)
                        [~,o1] = RobustAverage1D(V{g}(voxid,:)',option.robust_std);
                   else
                        o1 = zeros(size(V{g},2),1);
                   end
                   outliers = cat(1,outliers,o1);
                   Vl = cat(2,Vl,V{g}(voxid,:));
               end
%                [b,~,stats] = glmfit(X(~outliers,:),Vl(~outliers),'normal','constant','off');
%                [b,~,stats] = glmfit(X,Vl,'normal','constant','off');
               b = X(~outliers,:)\Vl(~outliers)';
               regression(voxid,:) = b;
%                regression_pval(voxid,:) = stats.p;
            end
        end
        
        nii.img = reshape(regression,sx,sy,sz,size(regression,2));
        
        if(option.smooth_std > 0)
            for dim4d=1:size(nii.img,4)
                V = nii.img(:,:,:,dim4);
                V = smooth3(V,'gaussian',option.smooth_std);
                nii.img(:,:,:,dim4) = V;
            end
        end        
        nii.hdr.dime.dim(1) = 4;
        nii.hdr.dime.dim(5) = size(nii.img,4);
        save_untouch_nii(nii,[templatesPath 'L_' num2str(l) '_RISHGLM_beta.nii.gz']);
        
        nii.img = reshape(regression_pval,sx,sy,sz,size(regression_pval,2));
        save_untouch_nii(nii,[templatesPath 'L_' num2str(l) '_RISHGLM_pval.nii.gz']);
    end

    % %% Difference Rish
    %
    % % For traveling heads,  subtract other site from reference site for each
    % % subject. Then take the mean.
    
    
    if option.travelingHeads
        % Not implemented
        
    else
        % For others, subtract the mean of other site from the mean of
        % reference site.
        for l=0:2:SHorder
            
            NGroups = length(SITES);
            
            refImg=load_untouch_nii([templatesPath  'Mean_' SITES{1}.name '_L'  num2str(l)  '.nii.gz']);
            movImg=load_untouch_nii([templatesPath  'Mean_' SITES{1}.name '_L'  num2str(l)  '.nii.gz']);
            
            rishglm = load_untouch_nii([templatesPath 'L_' num2str(l) '_RISHGLM_beta.nii.gz']);
            
            ref = option.reference_site;
            if(ref == -1)
                AVG_L = mean(rishglm.img(:,:,:,1:NGroups),4);
            else
                AVG_L = rishglm.img(:,:,:,ref);
            end

            for g=1:NGroups
    
                delta1 = rishglm.img(:,:,:,g) - AVG_L;
    
                scale1 = sqrt(AVG_L./(rishglm.img(:,:,:,g)+eps));
                scale1(scale1<0 | scale1 > 10) = 1;
                
                refImg.img = scale1.*single(maskRef.img);
                movImg.img = delta1.*single(maskRef.img);
                save_untouch_nii(refImg,[templatesPath  'Scale_L'  num2str(l) '_S' num2str(g) '.nii.gz']);
                save_untouch_nii(movImg,[templatesPath  'Delta_L'  num2str(l) '_S' num2str(g) '.nii.gz']);
                
            end
        end
               
        
    end
end    

end

function [medv,bV] = RobustAverage(Vol4D,std_tol)
    medv = zeros(size(Vol4D(:,:,:,1)));
    [sx,sy,sz,st] = size(Vol4D);
    UV = reshape(Vol4D,sx*sy*sz,st);
    mUV = max(UV,[],2);
    for x=1:size(UV,1)
       V = UV(x,:);
       if(mUV(x) == 0)
           continue
       end
       medianV = median(V);
       rSTD = 1.4826*mad(V,1);
       bV = abs(V-medianV)>=std_tol*rSTD;
       medv(x) = mean(V(bV==0));
    end
end

function [medv,bV] = RobustAverage1D(Signal1D,std_tol)
   medianV = median(Signal1D);
   rSTD = 1.4826*mad(Signal1D,1);
   bV = abs(Signal1D-medianV)>=std_tol*rSTD;
   medv = mean(Signal1D(bV==0));
end
