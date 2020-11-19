clear all
close all
clc

addpath(genpath(fullfile('..','Utils')))
Subj                        = '';

sensor_positions_file       = fullfile('..','Data',Subj,'Photogrammetry','sensor_positions.mat');
photogrammetry_mesh_file    = ''; % mesh of the subject derived from the 
% photogrammetry-based method
rescale_factor_file         = ''; % scaling factor 
scalp_surface_file          = ''; % mesh of the subject derived from T1 
% image (output of the function CreateVolumetricMesh.m)
mask_file                   = ''; % segmentation of the T1 image (nifti 
% file output of the function CreateMask.m)

load(sensor_positions_file)
load(rescale_factor_file)
load(scalp_surface_file)


ptCloud = pcread(photogrammetry_mesh_file);
loc     = double(ptCloud.Location)*ratio; 
scalp   = ScalpSurfaceMesh.node;
moving  = pointCloud(loc);
fixed   = pointCloud(scalp);

% Downsample photogrammetry-based mesh
range_pos       = max(moving.Location)-min(moving.Location);
grid_step       = double(mean(range_pos/30));
moving_low      = pcdownsample(moving,'GridAverage',grid_step);
node_moving_low = moving_low.Location;
face_moving_low = MyCrustOpen(node_moving_low);

% Downsample T1-based mesh
range_pos       = max(fixed.Location)-min(fixed.Location);
grid_step       = double(mean(range_pos/30));
fixed_low       = pcdownsample(fixed,'GridAverage',grid_step);
node_fixed_low  = fixed_low.Location;
face_fixed_low  = MyCrustOpen(node_fixed_low);

subj_p.node     = node_moving_low;
subj_p.face     = face_moving_low;
subj_mr.node    = node_fixed_low;
subj_mr.face    = face_fixed_low; 

% Select three points around the two ears and the nose in both surfaces.
% These points will be used to estimate a rigid transformation that 
% roughly aligned the two meshes in the same space and direction
heads           = {'subj_p';'subj_mr'};
for p = 1:2
    clear infro_struct
    clear dcm_obj
    fig = figure;
    eval(['node = ',heads{p},'.node;'])
    eval(['face = ',heads{p},'.face;'])
    plotmesh(node,face)
    disp('Move/rotate/zoom the figure, then press "Return"')
    pause
    datacursormode on
    dcm_obj = datacursormode(fig);
    for i=1:3
        disp('Click on the point, then press "Return"')
        pause 
        info_struct(i) = getCursorInfo(dcm_obj);
        if isfield(info_struct, 'Position')
          disp('Clicked positioin is')
          disp(info_struct(i).Position)
        end
    end
    datacursormode off
    eval(['p',num2str(p),' = [info_struct(:).Position];'])
end
p1 = reshape(p1,3,3)';
p2 = reshape(p2,3,3)';

figure,plot3(node_moving_low(:,1),node_moving_low(:,2),node_moving_low(:,3),'.k')
hold on,plot3(p1(:,1),p1(:,2),p1(:,3),'*r','MarkerSize',10)
figure,plot3(node_fixed_low(:,1),node_fixed_low(:,2),node_fixed_low(:,3),'.k')
hold on,plot3(p2(:,1),p2(:,2),p2(:,3),'*r','MarkerSize',10)

% Estimate and apply rigid transformation
[ret_R, ret_t] = rigid_transform_3D(p1,p2);
A2             = (ret_R*subj_p.node') + repmat(ret_t, 1, size(subj_p.node,1));
A2             = A2';

figure
subplot(2,1,1),plot3(A2(:,1),A2(:,2),A2(:,3),'.k')
subplot(2,1,2),plot3(node_fixed_low(:,1),node_fixed_low(:,2),node_fixed_low(:,3),'.k')

% Estimate and apply an affine transformation
figure
node_moving_low = A2;
clear opt Transform
opt.method      = 'affine'; % use affine registration
opt.viz         = 1;          % show every iteration
opt.normalize   = 1;    % normalize to unit variance and zero mean before registering (default)
opt.scale       = 1;        % estimate global scaling too (default)
opt.rot         = 1;          % estimate strictly rotational matrix (default)
opt.corresp     = 0;      % do not compute the correspondence vector at the end of registration (default)
opt.tol         = 1e-8;       % tolerance
opt.max_it      = 1000;
opt.outliers    = 0.8;
% registering Y to X
[Transform, Correspondence]=cpd_register(node_fixed_low,node_moving_low,opt);


for i = 1:length(s10_10_pos)
    pos_10_10(i,:) = s10_10_pos(i).pos;
end
for i = 1:length(s10_20_pos)
    pos_10_20(i,:) = s10_20_pos(i).pos;
end

pos     = [pos_10_10; pos_10_20; lm_pos(5).pos]*ratio;
pos_r   = (ret_R*pos') + repmat(ret_t, 1, size(pos,1));
pos_r   = pos_r';
pos_af  = (((Transform.R*pos_r')'+repmat(Transform.t',size(pos_r,1),1))).*Transform.s;


figure,
plot3(node_fixed_low(:,1),node_fixed_low(:,2),node_fixed_low(:,3),'.k')
hold on, plot3(pos_af(:,1),pos_af(:,2),pos_af(:,3),'*b')


for i = 1:size(pos_af,1)
    d               = sqrt(sum((pos_af(i,:)-ScalpSurfaceMesh.node).^2,2));
    [v,idx]         = min(d);
    new_pos_af(i,:) = ScalpSurfaceMesh.node(idx,:);
end

figure,
plotmesh(ScalpSurfaceMesh.node,ScalpSurfaceMesh.face,'FaceColor',[0.5 0.5 0.5])
hold on, plot3(new_pos_af(:,1),new_pos_af(:,2),new_pos_af(:,3),'*b','MarkerSize',20)

% Save the warped positions
pos         = round(new_pos_af);
save_file   = fullfile('..','Data',Subj,'T1','sensor_positions2T1.mat');
save(save_file,'pos')

% Create nifti image with the warped positions
mask_nii = load_untouch_nii(mask_file);
mask_img = mask_nii.img;
img      = zeros(size(mask_img));
for i = 1:size(pos,1)
    img(pos(i,1),pos(i,2),pos(i,3)) = i;
end
mask_nii.img = img;
save_file   = fullfile('..','Data',Subj,'T1','sensor_positions2T1.nii');
save_untouch_nii(mask_nii,save_file)



