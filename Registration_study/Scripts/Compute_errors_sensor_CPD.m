clear all
clc

addpath(genpath('../Utils/'))
load('../Data/MNI/ScalpSurfaceMesh.mat')
load('../Data/MNI/GMSurfaceMeshMNI.mat')

subj    = dir('../Data/Subj*');
transf  = 'CPD_non_linear'; %CPD_affine

for ss = 1:length(subj)
    load(fullfile(subj(ss).folder,subj(ss).name,transf,'Sensor_positions',[subj(ss).name,'2MNI.mat']))
    dataimg     = load_untouch_nii(fullfile(subj(ss).folder,subj(ss).name,'T1','SensorPos2MNI.nii'));  
    aimg        = double(dataimg.img);

    for i = 1:size(pos,1)
        xx              = find(aimg == i);
        [rr,cc,vv]      = ind2sub(size(aimg),xx);   
        temp            = mean([rr cc vv],1) - [96 132 78]; 
        [vv,ii]         = min(sqrt(sum(((temp - ScalpSurfaceMesh.node).^2),2)));
        true_pos (i,:)  = ScalpSurfaceMesh.node(ii,:);
    end
       
    euclidean_dist(:,ss) = sqrt(sum(((true_pos - pos).^2),2));

    for i = 1:size(pos,1)
        p1          = true_pos(i,:);
        p2          = pos(i,:);
        temp_p3     = [p1 + p2]/2;
        [vv,idx]    = min(sqrt(sum(((temp_p3 - ScalpSurfaceMesh.node).^2),2)));
        p3          = ScalpSurfaceMesh.node(idx,:);
        if norm(p1 - p3) == 0 | norm(p2 - p3) == 0
            geodesic_dist(i,ss) = euclidean_dist(i,ss);
        else
            [curve_pts_p_p1 len_p_p1] = curve_gen(p1,p2,p3,ScalpSurfaceMesh.node);
            if len_p_p1 < euclidean_dist(i,ss)
                geodesic_dist(i,ss) = euclidean_dist(i,ss);
            else
                geodesic_dist(i,ss) = len_p_p1;
            end
        end
    end
end

save_file = ['../Results/Sensor_pos_',transf,'.mat'];
save(save_file,'euclidean_dist','geodesic_dist')