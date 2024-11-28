% An L2 regularized NNLS function
function [Deconv_1,Deconv_1_clean] = RegularizedLS(Dictionary_1,OneBigVoxelFull_1,options,lambda,x0)
    if(nargin < 3)
        options = optimset('TolX',1e-12,'TolFun',1e-12);
    end

    % base_deconv = lsqnonneg(Dictionary_1,OneBigVoxelFull_1,options);
    % Amp = sum(base_deconv);
    Amp = 1;

    C = Amp*lambda*eye(size(Dictionary_1,2),size(Dictionary_1,2));
    C = C(1:size(Dictionary_1,2),1:size(Dictionary_1,2));
    Dictionary_1_padded = [Dictionary_1; C];
    OneBigVoxelFull_1_tmp = zeros(size(Dictionary_1_padded,1),1);
    OneBigVoxelFull_1_tmp(1:length(OneBigVoxelFull_1)) = OneBigVoxelFull_1;
    if(exist('x0','var') > 0)
        OneBigVoxelFull_1_tmp(length(OneBigVoxelFull_1)+1:end) = x0;
    end

    Deconv_1 = Dictionary_1_padded\double(OneBigVoxelFull_1_tmp);

    if(nargout > 1)
        Deconv_1_clean = zeros(size(Deconv_1));
        idx = Deconv_1>0;
        Deconv_1_clean(idx) = lsqnonneg(Dictionary_1(:,idx),OneBigVoxelFull_1,options);
    end

end