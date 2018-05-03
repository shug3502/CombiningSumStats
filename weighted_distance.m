function d=weighted_distance(y1,y2,weights)
%now assume y1 is Nxnum_paramsxrepeats
%weighted euclidean distance. Could probably use dist or an inbuilt fn
d = sum(sum((repmat(weights,size(y1,1),1,size(y1,3)).*bsxfun(@minus,y1,y2)).^2,2),3);
