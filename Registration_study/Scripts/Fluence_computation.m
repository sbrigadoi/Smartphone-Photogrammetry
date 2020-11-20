% This function allows the computation of the photon fluence on the MNI 
% template (case 1) and on individual subject anatomy (case 2).
% Case 1: load scalpSurfaceMesh.mat, TissueMask.mat (both located inside 
% the directory called Fluence_MNI) and the individual sensor positions 
% warped to the MNI. Sensor positions must be converted to the 
% scalp/tissue_mask space.
% When individual sensor positions are warped to the MNI template using 
% Affine_landmarks or CPD-based transformations use the following code:
% pos = round(sensor_positions +[96 132 78]); 
% When individual sensor positions are warped to the MNI template using 
% ANTs-based transformation use the following code:
% ras2lps  = diag([-1 -1 1 1]);
% mni_node = scalpSurfaceMesh.node-[96 132 178];
% pos  = xform_apply(dati,inv(ras2lps));
% for ii = 1:size(pos,1)
%    [vv,iii]    = min(sqrt(sum(((pos(ii,:) - mni_node).^2),2)));
%    pos(ii,:) = mni_node.node(iii,:);
% end
% pos=round(pos+[96 132 78]);
% Case 2: load the individual tissue mask (output of the function CreateMask.m),
% the individual scalp mesh (output of the function CreateVolumetricMesh.m), 
% and the sensor positions warped to the individual structural image 
%(output of the function Align_photogrammetry_mesh_2_T1_mesh).
mask_file   = '';
sensor_file = '';
scalp_file  = '';


wavelengths = 800;

load('Tissue_Properties_Model.mat');

anisotropy = [0.89 0.84 0.89 0.89 0.89];
refractive = repmat(1.3,1,5);
[~,wav_ind] = min(abs(lambda - wavelengths));
mua_vec = mua_model(wav_ind,:); 
mus_vec = mus_model(wav_ind,:)./(1-anisotropy); %reduced scattering coefficient

% Set MCX parameters
cfg.nphoton = 1e9;
%cfg.prop = [0 0 1 1; mua_vec' 0.02 9.0 0.89 1.37; 0.08 40.9 0.84 1.37; 0.004 0.009 0.89 1.37; ...
%    0.019 7.8 0.89 1.37; 0.019 7.8 0.89 1.37; 0 0 1 1];
cfg.prop = [0 0 1 1; mua_vec' mus_vec' anisotropy' refractive'; 0 0 1 1]; %air as background
% an N by 4 array, each row specifies [mua, mus, g, n] in order.
% the first row corresponds to medium type 0 which is
% typically [0 0 1 1]. The second row is type 1, and so on.
cfg.tstart= 0;
cfg.tstep = 5e-9;
cfg.tend = 5e-9;

%== MC simulation settings ==
cfg.seed=hex2dec('623F9A9E');
cfg.unitinmm = 1;
cfg.isnormalized = 1; % normalize the output fluence to unitary source

%== GPU settings ==
cfg.autopilot = 1;
cfg.gpuid=1;

%== Source-detector parameters ==
cfg.srctype = 'disk';
cfg.srcparam1 = 2; 
cfg.issrcfrom0 = 0;

%== Output control ==
cfg.outputtype = 'fluence'; % fluence integrated over each time gate,

%== Debug ==
cfg.debuglevel = 'P'; 


for iS = 1:nSubj
    
    %fprintf(['Working on Subject ' subj(iS).name '... \n'])
    pathnameSubj = fullfile(pathname,subj(iS).name);
	
    % Load mask
    load(fullfile(pathnameSubj,mask_file));  
	
    % Load sensor positions
    load(fullfile(pathnameSubj,sensor_file)); 

    load(fullfile(pathnameSubj,scalp_file)); 
    % Compute the center of the head   
    center_of_head  = mean(scalp,1);
    
    % Find the vector between the source/detector and the center of the head
    % and normalize it
    fprintf('Computing the normal vectors... \n')
    vector_src = repmat(center_of_head,size(pos,1),1) - pos;
    vector_norm_src = zeros(size(vector_src));
    for i_r = 1:size(vector_src,1)
        vector_norm_src(i_r,:) = vector_src(i_r,:)./norm(vector_src(i_r,:));
    end
    
    cfg.vol = uint8(mask);
   
    for pp = 1:size(pos,1)
        cfg.srcpos = pos(pp,:);
        cfg.srcdir = vector_norm_src(pp,:);
        fluenceSrc(pp) = mcxlab(cfg);
    end
    save(fullfile(pathnameSubj,'fluenceSrc.mat'),'fluenceSrc','-v7.3');
end