clear all
clc

addpath(genpath('../Utils/'))
load('../Other/GMSurfaceMesh.mat')
load('../Other/phantom_positions.mat')

[num,txt,raw]   = xlsread('../Other/name_positions.xlsx');
label_ref       = {txt{2:289,1}}; 

subj            = dir('../Data/Subj*');
transf          = {'Affine_landmarks','CPD_affine','CPD_non_linear','ANTs_affine','ANTs_non_linear'};
count           = [];

for ss = 1:length(subj) 
    pathfile = fullfile(subj(ss).folder,subj(ss).name,'T1','Act_loc2MNI.mat');
    load(pathfile)
    true_pos = loc;
    
    for tt = 1:length(transf)
        pathfile                    = fullfile(subj(ss).folder,subj(ss).name,transf{tt},'Fluence','Act_loc2MNI.mat');
        load(pathfile)
        new_pos                     = loc;  
        
        euclidean_dist_all(tt,:,ss) = sqrt(sum(((true_pos - new_pos).^2),2));

        for i = 1:size(new_pos,1)
            p1          = true_pos(i,:);
            p2          = new_pos(i,:);
            temp_p3     = [p1 + p2]/2;
            [vv,idx]    = min(sqrt(sum(((temp_p3 - GMSurfaceMesh.node).^2),2)));
            p3          = GMSurfaceMesh.node(idx,:);
            if norm(p1 - p3) == 0 | norm(p2 - p3) == 0
                geodesic_dist_all(tt,i,ss) = euclidean_dist_all(tt,i,ss);
            else
                [curve_pts_p_p1 len_p_p1] = curve_gen(p1,p2,p3,GMSurfaceMesh.node);
                if len_p_p1 < euclidean_dist_all(tt,i,ss)
                    geodesic_dist_all(tt,i,ss) = euclidean_dist_all(tt,i,ss);
                else
                    geodesic_dist_all(tt,i,ss) = len_p_p1;
                end
            end
        end
        clear loc new_pos
    end
    clear true_pos loc
end



for tt = 1:length(transf) 
    geodesic_dist   = squeeze(geodesic_dist_all(tt,:,:));
    euclidean_dist  = squeeze(euclidean_dist_all(tt,:,:));
    save_file       = ['../Results/Fluence_act_pos_',transf{tt},'.mat'];
    save(save_file,'euclidean_dist','geodesic_dist')
end
        
        

    
    
    
    
    
    
    
    
    
    
    
    

        
        
        

