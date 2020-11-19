clear all
close all
clc


addpath(genpath(fullfile('..','Utils')))

subj           = 'Subj1';

load(fullfile('..','Other','ScalpSurfaceMesh.mat'))
load(fullfile('..','Other','phantom_positions.mat'))
load(fullfile('..','Other','HeadVolumeMesh.mat'))
load(fullfile('..','Other','GMSurfaceMeshMNI.mat'))
load(fullfile('..','Data',subj,'Photogrammetry','sensor_positions.mat'));
photogrammetry_mesh_file    = ''; % mesh of the subject derived from the 
% photogrammetry-based method

[num,txt,raw]  = xlsread(fullfile('..','Other','name_positions.xlsx'));
label_ref      = {txt{2:289,1}}; 

subj_level     = -17; % level used to cut the individual mesh with an axial plane 
% under the nasion. This parameter must be manually set 
mni_level      = -55;

% Downsample the MNI surface mesh
fixed          = pointCloud(ScalpSurfaceMesh.node);
range_pos      = max(fixed.Location) - min(fixed.Location);
grid_step      = double(mean(range_pos/50));
fixed_low      = pcdownsample(fixed,'GridAverage',grid_step);
node_fixed_low = double(fixed_low.Location);
node_fixed_low(find(node_fixed_low(:,3)<mni_level),:) = []; % Cut the MNI surface
% with an axial plane under the nasion

lm_label = {lm_pos(:).name};
for i = 1:length(lm_label), 
    idx_lm(i)    = find(strcmpi(lm_label{i},label_ref)); 
end
for i = 1:length(s10_10_pos)
    pos_10_10(i,:) = s10_10_pos(i).pos;
end
for i = 1:length(s10_20_pos)
    pos_10_20(i,:) = s10_20_pos(i).pos;
end        
pos         = [pos_10_10;pos_10_20;lm_pos(5).pos];
marker      = [lm_pos(:).pos];
marker      = reshape(marker,3,5)';
p1a         = marker;
p1a(find(strcmp({lm_pos(:).name},'inion')),:) = []; % We do not need inion position
p2a         = surface_points(idx_lm,:);
idx_inion   = find(strcmpi(label_ref,'inion'));
p2a(find(idx_lm == idx_inion),:)              = [];
[ret_R, ret_t]  = rigid_transform_3D(p1a,p2a); % rigid transformation estimated using four easily detectable points 
% on a mesh, i.e., nasion, Cz, left and right preauricular points

pt         = pcread(photogrammetry_mesh_file);
moving     = pt.Location;
moving     = (ret_R*moving') + repmat(ret_t, 1, size(moving,1));
moving     = moving';
pos        = (ret_R*pos') + repmat(ret_t, 1, size(pos,1));
pos        = pos';
range_pos  = max(moving)-min(moving);
grid_step  = double(mean(range_pos/50));
moving_low = pcdownsample(pointCloud(moving),'GridAverage',grid_step);
node_moving_low = double(moving_low.Location);
node_moving_low(find(node_moving_low(:,3) < subj_level),:) = [];
face_moving_low      = MyCrustOpen(node_moving_low);
[conn,connnum,count] = meshconn(face_moving_low,size(node_moving_low,1));
node_moving_low      = smoothsurf(node_moving_low,[],conn,1000,0.5,'lowpass');
    
% For each sensor position find the nearest point on the surface 
for i = 1:size(pos,1)
    [aa,bb]             = min(sqrt(sum(((pos(i,:)-node_moving_low).^2),2)));
    idx_pos_on_surf(i)  = bb;
end

clear opt Transform
opt.method      = 'nonrigid'; % use rigid registration
opt.viz         = 1; % show every iteration
opt.beta        = 4;
opt.lambda      = 5;
opt.normalize   = 1;    % normalize to unit variance and zero mean before registering (default)
opt.scale       = 1;        % estimate global scaling too (default)
opt.rot         = 1;          % estimate strictly rotational matrix (default)
opt.corresp     = 0;      % do not compute the correspondence vector at the end of registration (default)
opt.tol         = 1e-8;       % tolerance
opt.max_it      = 2000;
opt.outliers    = 0.4;
% registering Y to X
[Transform, Correspondence] = cpd_register(node_fixed_low,node_moving_low,opt);
    
pos = Transform.Y(idx_pos_on_surf,:);

figure,
plotmesh(ScalpSurfaceMesh.node,'.k'), hold on,
plotmesh(pos,'*r','MarkerSize',20)

id       = subj(end);
savefile = fullfile('..','Data',subj,'CPD_non_linear',['Subj',id,'2MNI.mat']);
save(savefile,'pos')


        
    
    
    

    
    
    
    
 

