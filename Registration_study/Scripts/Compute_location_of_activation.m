clear all
close all
clc


addpath(genpath('../Utils/'))
load('../Other/GMSurfaceMesh.mat')

subj            = dir('../Data/Subj*');
transf          = {'Affine_landmarks','CPD_affine','CPD_non_linear','ANTs_affine','ANTs_non_linear'};
 
thresh          = 0.8;

for ss = 1:length(subj)
    
    for tt = 1:length(transf)
        pathfile = fullfile(subj(ss).folder,subj(ss).name,transf{tt},'Fluence');
        load (fullfile(pathfile,'Fluence2MNI.mat'))
   
        for ch = 1:size(fluence,2)
            fl          = fluence(:,ch);
            val         = max(fl);
            idx_max     = find(fl>thresh*val);
            pos_max     = GMSurfaceMesh.node(idx_max,:);
            val_max     = fl(idx_max);
            w           = val_max/val;      
            loc(ch,:)   = sum(pos_max.*w,1)/sum(w);
            clear fl val idx_max pos_max val_max w
        end
        
        savefile = fullfile(subj(ss).folder,subj(ss).name,transf{tt},'Fluence','Act_loc2MNI.mat');
        save(savefile,'loc')
        clear loc
    end
end