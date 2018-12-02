function R = my_reshape(M)
%% reshape not quite how I need it
%%%%%%%%
sz = size(M);
if length(sz)==3
    R = zeros(sz(1)*sz(3),sz(2));
    for k=1:sz(3)
        R((k-1)*sz(1)+(1:sz(1)),:) = M(:,:,k);
    end
else 
    R = M;
end
