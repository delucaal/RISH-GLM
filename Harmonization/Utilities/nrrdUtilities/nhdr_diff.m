function [S u b voxel baseline] = nhdr_diff(fn)
  paths
  
  dwi = loadNrrdStructure(fn);
  S = dwi.data;
  u = dwi.gradients;
  M = dwi.measurementframe;
  b = dwi.bvalue;
  voxel = abs(diag(dwi.spacedirections));
  
  ndirs = size(u,1);
  
  % first find the number of baseline data
  ct = 1;
  while norm(u(ct,:)) == 0
    ct = ct + 1;
  end
  
  ndir = ndirs - ct + 1;
  
  % separate out baseline
  s0 = S(:,:,:,1:ct-1);
  S  = S(:,:,:,ct:end); % drop null-gradient slices
  
  % actual gradient directions used
  u = u(ct:end,:);
  u = [u; -u]; % antipodal
  
  % divide off baseline
  s0 = mean(s0, 4); if nargout == 4, baseline = s0; end
  
  s0(s0==0) = 1; 
  S = single(S) ./ s0(:,:,:,ones(1,ndir));
  
  %S = single(S) ./ repmat(s0,[1 1 1 ndir]);
  
  
end
