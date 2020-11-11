% This function applies an affine transformation matrix T to each point and
% outputs the transformed array of points

function ptsTransf = applyTransf(pts,T)

n = size(pts, 1);
pts=pts';

if(all(size(T) == [4 4]) || all(size(T) == [3 4]))
    A = T(1:3, 1:3);
    b = T(1:3,4);
elseif(all(size(T) == [3 3]) || all(size(T) == [2 3]))
    A = T(1:2,1:2);
    b = T(1:2,3);
end

ptsTransf = A*pts;
for iP = 1:n
    ptsTransf(:,iP) = ptsTransf(:,iP)+b;
end
ptsTransf = ptsTransf';
