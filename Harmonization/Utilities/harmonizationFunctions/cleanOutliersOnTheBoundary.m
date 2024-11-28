function ScaleLn =  cleanOutliersOnTheBoundary(maskring,scalemap, threshold, flag)

ring = single(maskring).*scalemap;
% nm  = ring>0;
% idx = kmeans(ring(nm), 2);
% 
% muCl1 = mean(ring(idx==1));
% muCl2 = mean(ring(idx==2));
% 
% if muCl1>muCl2
%    threshold = min(ring(idx==1));
% else
%    threshold = min(ring(idx==2));
% end
% 
cleanedScaleMap = scalemap;

if threshold>10
    cleanedScaleMap(scalemap>prctile(scalemap(maskring>0), 95))=1;%prctile(scalemap(maskring>0), 95);
    threshold = prctile(scalemap(maskring>0), 95);
end

if flag % boundary
    cleanedScaleMap(ring>=threshold)=eps;
end
ind=find(ring>=threshold);
ScaleLn = medfilt_roi(cleanedScaleMap, ind);
end