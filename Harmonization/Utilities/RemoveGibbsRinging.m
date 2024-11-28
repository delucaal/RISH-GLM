%Works in matlab 2013a/2015b/2017a
function  deGibbs = RemoveGibbsRinging(S)

iter = 1;%number of times to repeat unringing software

addpath ./Utilities/unring/matlab


deGibbs = unring(S);

%if we need to repeat the unringing for severe Gibbs ringing
if iter >= 2
    for i=2:iter
      deGibbs = unring(deGibbs);
    end
end


end

