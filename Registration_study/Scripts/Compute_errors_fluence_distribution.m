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
    pathfile        = fullfile(subj(ss).folder,subj(ss).name,'T1','Fluence2MNI.mat');
    load(pathfile)
    true_fl         = fluence;
	
	for tt = 1:length(transf)
        pathfile = fullfile(subj(ss).folder,subj(ss).name,transf{tt},'Fluence','Fluence2MNI.mat');
        load(pathfile)
        est_fl              = fluence;
			
		for ch = 1:size(est_fl,2)
            f1   = true_fl(:,ch);
            val  = max(f1);
            val1 = zeros(size(f1));
            idx_max = find(f1 > 0.80*val);
            val1(idx_max)  = 1;
            true_fl_points = GMSurfaceMesh.node(idx_max,:);

            f1   = est_fl(:,ch);
            val  = max(f1);
            val2 = zeros(size(f1));
            idx_max = find(f1 > 0.80*val);
            val2(idx_max)  = 1;
            est_fl_points  = GMSurfaceMesh.node(idx_max,:); 
        
            dice1(ch) = dice(val1,val2);  
        end
     	
        dice_all(tt,ss,:)   = dice1;
   end
end

for tt = 1:length(transf) 
    dice     = squeeze(dice_all(tt,:,:));
    save_file       = ['../Results/Fluence_distrib_dice_',transf{tt},'.mat'];
    save(save_file,'dice')
end






