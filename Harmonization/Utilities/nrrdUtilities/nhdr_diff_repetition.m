function nhdr_diff_repetition(fn,N)

 %output
 %input, fn - file name
 %N - spacing between each repetitiion
 
  paths;
  
  dwi = loadNrrdStructure(fn);
  S = dwi.data;
  %u = dwi.gradients;
  if isfield(dwi, 'measurementframe')
    M = dwi.measurementframe;
  else
    M = eye(3);
  end
  b = dwi.bvalue;
  order = [1 2 3 4];
  %make adjustment of the spacedirections and measurement frame
  
  sd = dwi.spacedirections;
  voxel = [norm(sd(:,1));norm(sd(:,2));norm(sd(:,3))]';
  R = sd./repmat(voxel,[3 1]);
  
  %u = M*inv(R)*dwi.gradients';
  %u = u';
 
  u = inv(R)*M*dwi.gradients';u=u';
  for j = 1:size(u,1)
      if norm(u(j,:)) > eps
       u(j,:) = u(j,:)/norm(u(j,:));
      end
  end
  %Sm is the space directions needed to transform the fiber tracts in the
  %slicer space
  if strcmp(dwi.space,'LPS') || strcmpi(dwi.space,'left-posterior-superior')
    RAS = [-1 0 0;0 -1 0;0 0 1];
    Sm = RAS * dwi.spacedirections;
  else
    Sm = dwi.spacedirections;
  end
  %Sm = inv(R)*M;
  
  ndirs = size(u,1);
  
  %find the storage format
  switch size(u,1)
      case size(S,1)
         order = [2 3 4 1];
      case size(S,2)
         order = [1 4 3 2];
      case size(S,3)
          order = [1 2 4 3];
  end
  S = permute(S,order);
  sz = size(S);
  Sm(4,:) = -0.5 * Sm * sz(1:3)'; %the space origin in Slicer
  
  
  % first find the number of baseline data
  idb=[];
  id=[];
  for j =1:ndirs
      %the baseline images
      if norm(u(j,:)) == 0
          idb = [idb;j];
      else
          id = [id;j];
      end
  end
  
  if isempty(idb)
      fprintf('Data does not have a B0 image. \n');
      return;
  end
  ndir = length(id);
  
 
  % separate out baseline
  s0 = S(:,:,:,idb);
  S  = S(:,:,:,id); % drop null-gradient slices
  
  % actual gradient directions used
  u = u(id,:);
  u = [u; -u]; % antipodal
  
  % divide off baseline
  s0 = mean(s0, 4); if nargout == 4, baseline = s0; end
  
  s0(s0==0) = 1; 
  S = single(S) ./ s0(:,:,:,ones(1,ndir));
  
  
  
end
