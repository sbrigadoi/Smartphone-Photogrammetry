clear all
clc

addpath(genpath('../Utils/'))
load('../Data/MNI/ScalpSurfaceMesh.mat')
load('../Data/MNI/GMSurfaceMeshMNI.mat')
load('../Other/phantom_positions.mat')

[num,txt,raw]   = xlsread('../Other/name_positions.xlsx');
label_ref       = {txt{2:289,1}}; 

subj            = dir('../Data/Subj*');
transf          = 'Affine_landmarks'; 

for ss = 1:length(subj)
    
    clear idx_10_20 idx_10_10 idx_lm
    
    dataimg     = load_untouch_nii(fullfile(subj(ss).folder,subj(ss).name,'T1','SensorPos2MNI.nii'));  
    aimg        = double(dataimg.img);
    
    load(fullfile(subj(ss).folder,subj(ss).name,'Photogrammetry','sensor_positions.mat'))
    
    for i = 1:length(lm_pos)      
        idx_lm(i)      = find(strcmpi(lm_pos(i).name,label_ref)); 
        lm(i,:)        = lm_pos(i).pos;
    end
    for i = 1:length(s10_10_pos)      
        idx_10_10(i)   = find(strcmpi(s10_10_pos(i).name,label_ref)); 
        pos_10_10(i,:) = s10_10_pos(i).pos;
    end
    for i = 1:length(s10_20_pos)      
        idx_10_20(i)   = find(strcmpi(s10_20_pos(i).name,label_ref)); 
        pos_10_20(i,:) = s10_20_pos(i).pos;
    end

    pos     = [pos_10_10; pos_10_20; lm(5,:)];
    p1      = lm;
    p2      = surface_points(idx_lm,:);
    T       = gen_xform_from_pts(p1,p2);
    new_pos = xform_apply(pos,T);
    
    for i = 1:size(pos,1)
        xx              = find(aimg == i);
        [rr,cc,vv]      = ind2sub(size(aimg),xx);   
        temp            = mean([rr cc vv],1) - [96 132 78]; 
        [vv,ii]         = min(sqrt(sum(((temp - ScalpSurfaceMesh.node).^2),2)));
        true_pos (i,:)  = ScalpSurfaceMesh.node(ii,:);
    end
        
    for i = 1:size(pos,1)
        [vv,ii]         = min(sqrt(sum(((new_pos(i,:) - ScalpSurfaceMesh.node).^2),2)));
        new_pos(i,:)    = ScalpSurfaceMesh.node(ii,:);
    end
    
    euclidean_dist(:,ss) = sqrt(sum(((true_pos - new_pos).^2),2));
    
    for i = 1:size(pos,1)
        p1 = true_pos(i,:);
        p2 = new_pos(i,:);
        temp_p3   = [p1 + p2]/2;
        [vv,idx]  = min(sqrt(sum(((temp_p3 - ScalpSurfaceMesh.node).^2),2)));
        p3        = ScalpSurfaceMesh.node(idx,:);
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