% This function computes the affine transformation matrix from a set of
% point to another set of points

function T = affineTransf(p1,p2)

[p,~] = size(p1);
[~,n] = size(p2);

T = eye(n+1);
A = [p1, ones(p,1)];
for iP=1:n
    x(:,iP) = pinv(A)*p2(:,iP);
    T(iP,:) = x(:,iP)'; 
end