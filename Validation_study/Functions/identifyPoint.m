% This function identifies the detected points belonging to the same cluster
% and output the cluster's center as the point coordinate
%
% INPUT: 
%          loc:         locations (n x 3) of all the identified points
%          distThresh:  distance within which the points are considered
%                       belonging to the same cluster
%
% OUPUT:
%          pos:         coordinates (n x 3) of the identified points

function pos = identifyPoint(loc,distThresh)

temp1 = nan(size(loc,1));
for i = 1:size(loc,1)
    current_pt  = loc(i,:);
    temp1(i,:)  = sqrt(sum((loc-current_pt).^2,2));
end
count = 1;
done  = [];
for i = 1:size(temp1,1)
    if sum(i==done) == 0
        current_dist    = temp1(i,:);
        nn              = find(current_dist<distThresh);
        new             = [];
        for j = 1:length(nn)
            new_dist    = temp1(nn(j),:);
            new_nn      = find(new_dist<distThresh);
            new         = [new new_nn];
        end
        new             = unique(new);
        nn              = unique([new nn]);
        cl(count).ngb   = nn;
        done            = [done nn];
        count           = count+1;
    end
end

% Find center of cluster
pos = nan(length(cl),3);
for i = 1:length(pos)
    pp              = cl(i).ngb;
    pp_pos          = loc(pp,:);
    pos(i,:)  = mean(pp_pos,1);
end