function [consistentDWI] = nrrdReformatAndNormalize(rawDWI)
 %  This function will review parameters from a dwi struct to ensure
 %  that all necessary information is provided, inject "resonable defaults" for missing
 %  values, and permute the data into a canonical organization for subsequent
 %  data processing.
 %
 % --The measurement frame is reset to voxel lattice orientation (and then
 %   set to identity)
 % --The gradient directions are normailzed
 % --The order is permuted to consistent space vs. gradient memory layout
 %
 % Author Hans J. Johnson

  consistentDWI = rawDWI;

  %XXXXXXXXXXXXX
  % Make an identity measurement frame by rotating gradients to be
  % interpreted in voxel lattice orientation.
  if ~isfield(consistentDWI, 'measurementframe')
    consistentDWI.measurementframe = eye(3);
  end


 % KW QUESTION: All the 'HACK' stuff in this file has never been resolved.
 % What is the correct thing to do here?
 %make adjustment of the spacedirections and measurement frame
 % HACK : I think that using the direction cosignes is wrong here.
 spaceDirectionMatrix = rawDWI.spacedirections;
 voxel = [norm(spaceDirectionMatrix(:,1));norm(spaceDirectionMatrix(:,2));norm(spaceDirectionMatrix(:,3))]';
 directionCosines = spaceDirectionMatrix./repmat(voxel,[3 1]);


  % Remove the measurement frame from the gradient direcitons.  This makes
  % gradients relative to the voxel lattice.
  consistentDWI.gradientdirections = ( inv(directionCosines)*consistentDWI.measurementframe*consistentDWI.gradientdirections' )';
  % Now force to an identity measurement frame
  %  HACK: I think the above computation with the directionCosigns is wrong,
  %  HACK: I also think that this line should be added back in
  %consistentDWI.measurementframe = eye(3);
  %XXXXXXXXXXXXX

  % Renormalize to ensure unit length gradients (force all bValues to be
  % the same!!
  for j = 1:size(consistentDWI.gradientdirections,1)
      gnorm = norm(consistentDWI.gradientdirections(j,:));
      if gnorm > eps
       consistentDWI.gradientdirections(j,:) = consistentDWI.gradientdirections(j,:)/gnorm;
       % Should issue a warning here, because it indicates new b0 values
      end
  end

  %XXXXXXXXXXXXX
  % Permute the order of data to be in a cononical format
  order = [1 2 3 4];
  %find the storage format
  switch size(consistentDWI.gradientdirections,1)
      case size(consistentDWI.data,1)
         order = [2 3 4 1];
      case size(consistentDWI.data,2)
         order = [1 4 3 2];
      case size(consistentDWI.data,3)
          order = [1 2 4 3];
  end
  consistentDWI.data = permute(consistentDWI.data,order);
  consistentDWI.centerings = consistentDWI.centerings(order);
  consistentDWI.kinds = consistentDWI.kinds(order);
  %XXXXXXXXXXXXX
end
