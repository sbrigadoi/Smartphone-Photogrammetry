clear all
close all
clc

addpath(genpath(fullfile('..','Utils')))

Subj     = '';
photogrammetry_mesh_file    = ''; %mesh of the subject derived from the photogrammetry based method. Mesh should be in .ply format
load  (fullfile('..','Data',Subj,'Photogrammetry','sensor_positions.mat'));
rescale_factor_file         = ''; % scaling factor
load(rescale_factor_file)

load (fullfile('..','Other','HeadVolumeMesh.mat'))
load (fullfile('..','Other','ScalpSurfaceMesh.mat'))
load (fullfile('..','Other','phantom_positions.mat'))
[num,txt,raw]   = xlsread(fullfile('..','Other','name_positions.xlsx'));

subj_level     = -17; % level used to cut the individual mesh with an axial plane 
% under the nasion. This parameter must be manually set 

% Estimate an affine transformation between landmarks to realign the
% subject mesh to the phantom mesh
label_ref       = {txt{2:289,1}};
lm_label        = {lm_pos(:).name};

for i = 1:length(lm_label)
    idx_lm(i) = find(strcmp(lm_label{i},label_ref));
    marker(i,:)      = lm_pos(i).pos;
end
for i = 1:length(s10_10_pos)
    pos_10_10(i,:) = s10_10_pos(i).pos;
end
for i = 1:length(s10_20_pos)
    pos_10_20(i,:) = s10_20_pos(i).pos;
end        

pos         = [pos_10_10;pos_10_20;lm_pos(5).pos];
p1a         = marker;
p1a(find(strcmp({lm_pos(:).name},'inion')),:) = []; % We do not need inion position
p2a         = surface_points(idx_lm,:);
idx_inion   = find(strcmpi(label_ref,'inion'));
p2a(find(idx_lm == idx_inion),:)              = [];
[ret_R, ret_t]  = rigid_transform_3D(p1a,p2a); % rigid transformation estimated using four easily detectable points 
% on a mesh, i.e., nasion, Cz, left and right preauricular points

ptCloud     = pcread(photogrammetry_mesh_file);
range_pos   = max(ptCloud.Location)-min(ptCloud.Location);
grid_step   = double(mean(range_pos/50));
temp        = pcdownsample(ptCloud,'GridAverage',grid_step);
moving      = double(temp.Location);

% Apply the rigid transformation to the subject surface
A2 = (ret_R*moving') + repmat(ret_t, 1, size(moving,1));
A2 = A2';
A3 = A2;
A3(find(A3(:,3)<subj_level),:) = []; % Cut the realigned surface below the nose
face                    = MyCrustOpen(A3); % Compute the triangulation and smooth the surface
[conn,connnum,count]    = meshconn(face,size(A3,1));
node                    = smoothsurf(A3,[],conn,1000,0.5,'lowpass'); 

loc1                    = double(node*ratio);
node2                   = loc1-mean(loc1)+mean(ScalpSurfaceMesh.node); % Align
% the centre of mass of the individual scalp surface with the centre of mass of the template surface
pos                     = (ret_R*pos') + repmat(ret_t, 1, size(pos,1));
pos                     = pos'*ratio;
pos                     = pos-mean(loc1)+mean(ScalpSurfaceMesh.node);
    
for ch = 1:size(pos,1)
    [xx,idx]    = min(sqrt(sum(((pos(ch,:)-node2).^2),2)));
    pos(ch,:)   = node2(idx,:);
end 

% Create the nifti files and csv file that will be used by the function
% Photogrammetry2MNI_ANTs_affine.m and Photogrammetry2MNI_ANTs_non_linear.m
step            = 1.5;
grid            = [-150:step:150];
[img5 v2smap]   = surf2vol(node2,face,grid,grid,grid);
imga_nii        = make_nii(img5,[step step step],[round(length(grid)/2) round(length(grid)/2) round(length(grid)/2)],16);
save_file       = fullfile('..','Data',subj,'Photogrammetry','Subj_surface_no_nose_ANTs.nii');
save_nii(imga_nii,save_file) 
    
label           = {'left';'nasion';'right';'cz'};
lll             = {lm_pos(:).name};
r               = 3;
area            = 0.5;
marker          = [lm_pos(:).pos];
marker          = reshape(marker,3,5)';
p1_rigid        = (ret_R*marker') + repmat(ret_t, 1, size(marker,1));
p1_rigid        = p1_rigid';  
marker          = p1_rigid*ratio-mean(loc1)+mean(ScalpSurfaceMesh.node);
    
figure
plot3(node2(:,1),node2(:,2),node2(:,3),'.k')
hold on
plot3(marker(:,1),marker(:,2),marker(:,3),'*c')
hold on
plot3(pos(:,1),pos(:,2),pos(:,3),'*r')
    
ras2lps         = diag([-1 -1 1  1]);
pos_lps         = xform_apply(pos,ras2lps);
save_file       = fullfile('..','Data',subj,'Photogrammetry','Sensor_positions_ANTs.csv');
csvwrite(save_file,pos_lps)
    

for i = 1:size(label)
    idx                 = find(strcmp(label{i},lll));
    pos                 = marker(idx,:);
    loc_sorted(i,:)     = pos;
    [cnode,cface]       = meshasphere(pos,r,area);
    [cnode,cface]       = meshcheckrepair(cnode,cface);
    eval(['[img',num2str(i),', v2smap] = surf2vol(cnode,cface,grid,grid,grid,''fill'',''1'');'])
end
img         = img1+2*img2+3*img3+4*img4;
imga_nii    = make_nii(img,[step step step],[round(length(grid)/2) round(length(grid)/2) round(length(grid)/2)],16);
save_file   = fullfile('..','Data',subj,'Photogrammetry','Subj_surface_no_nose_spheres_ANTs.nii');
save_nii(imga_nii,save_file)

    


    
    











