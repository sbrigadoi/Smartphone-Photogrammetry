%This function performs the mapping from a tetrahedral volume mesh to the
%associated GM mesh.  Output is in the form of a sparse transformation matrix with
%dimensions NxM where M is the number of tetrahedral mesh nodes and N is
%the number of GM surface nodes;

clear all
close all
clc


addpath(genpath('../Utils/'))
load('../Other/GMSurfaceMesh.mat')

subj            = dir('../Data/Subj*');
transf          = {'Affine_landmarks','CPD_affine','CPD_non_linear','ANTs_affine','ANTs_non_linear'};

for ss = 1:length(subj)
    
    for tt = 1:length(transf)
        pathfile = fullfile(subj(ss).folder,subj(ss).name,transf{tt},'Fluence');
        load (fullfile(pathfile,'fluenceSrc.mat'))
       
        for ch = 1:length(fluenceSrc)

            fluence = double(fluenceSrc(ch).data);
            for i = 1:size(GMSurfaceMesh.node,1)
                curr            = round(GMSurfaceMesh.node(i,:));
                lim2            = curr+1;
                lim1            = curr-1;
                fluence(i,ch)   = mean(mean(mean(fluence_lin(lim1(1,1):lim2(1,1), lim1(1,2):lim2(1,2), lim1(1,3):lim2(1,3)))));
            end
        end
        savefile = fullfile(subj(ss).folder,subj(ss).name,transf{tt},'Fluence','Fluence2MNI.mat');
        save(savefile,'fluence')
        clear fluence
    end 
end



    
  



 