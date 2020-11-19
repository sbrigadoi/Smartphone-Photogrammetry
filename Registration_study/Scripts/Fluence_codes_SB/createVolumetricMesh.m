function [mesh,quality] = createVolumetricMesh(mask,maxvol,voxel_size,opt)

% Function that creates a volumetric tetrahedral mesh from a segmented
% mask. It uses the iso2mesh toolbox. The output mesh is scaled in real
% proportion
%
% Input: mask             segmented mask; it can contain different tissue types, 
%                         which will be meshed keeping the tissue information 
%
%        maxvol           maximum volume that a tetrahedron can have
%
%        voxel_size       voxel size (from the original MRI from which the
%                         mask comes from). It is needed to scale the mesh
%                         to real dimension
%
%        opt              different option for the creation of the mesh,
%                         see the help of vol2surf for more information 
%
% Output: mesh            the created mesh. It is a struct with fields
%                         mesh.node (number of nodes x 4, where the first 3 columns represent 
%                         the nodes coordinates and the fourth a label),
%                         mesh.elem (n of elements x 5, where the first 4
%                         colums are the node index of the tetrahedron and
%                         the 5th one the tissue type) and mesh.face
%                         (number of faces x 4, where the first 3 columns
%                         are the index of the nodes forming the triangle
%                         and the 4th one the tissue type). Faces form all
%                         the surfaces of the mesh.
%                         
%         quality         a vector number of elem x 1, which contains
%                         values between 0 and 1. Values close to 1
%                         indicate a good mesh, values close to 0 a
%                         degenerated tetrahedron
%
% Written by S. Brigadoi, UCL, 03/2013
%
%

% In the paper we used:
%mask 			= MaskCut;
%voxel_size 	= 1;
%opt.maxnode 	= 1000000; 
%maxvol 		= 3;

dim = size(mask);

% Create tetrahedral mesh
tic
[node,elem,face] = vol2mesh(uint8(mask),1:dim(1),1:dim(2),1:dim(3),opt,maxvol,1,'cgalmesh');
toc
figure,plotmesh(node,face)
%%
% Reorient the element of the mesh
elem_or = meshreorient(node(:,1:3),elem(:,1:4)) ;
elem = [elem_or,elem(:,5)];

if length(voxel_size)>1
    voxel_size = repmat(voxel_size,size(node,1),1);
end

node(:,1:3) = node(:,1:3).*voxel_size;     % because the output of iso2mesh is in grid/voxel_size: in order 
                                           % to make the mesh having the same dimension
                                           % of the real MRI we need to multiply
                                           % for the voxel size;
mesh.node = node;
mesh.elem = elem;
mesh.face = face;

% Compute quality of the mesh
%quality = meshquality(node(:,1:3),elem(:,1:4));

% [conn,connnum,count]=meshconn(face,size(node,1));
% node=smoothsurf(node,[],conn,1000,0.5,'lowpass');
% figure, plotmesh(node,face)
%%
tiss = unique(mesh.elem(:,5)); %#ok<NODEF>

% Use histograms to compute how frequently a node is shared by elements in
% the same tissue and then assign to that node the tissue having the
% highest frequency
tmp = zeros(size(mesh.node,1),length(tiss));
for i_t = 1:length(tiss)
    tiss_nodes = mesh.elem(mesh.elem(:,5) == i_t,1:4);
    [n,~] = hist(tiss_nodes(:),1:size(mesh.node,1));
    tmp(unique(tiss_nodes(:)),i_t) = n(unique(tiss_nodes(:)));
end
[~,ind] = max(tmp,[],2);
mesh.node(:,4) = ind;

save mesh mesh
%%
type={'scalp'};
ScalpSurfaceMesh = surfaceMeshLD(mask,type,voxel_size);
scalp=ScalpSurfaceMesh.node;
save ScalpSurfaceMesh ScalpSurfaceMesh scalp
figure, plotmesh(ScalpSurfaceMesh.node,ScalpSurfaceMesh.face)
figure, plot3(ScalpSurfaceMesh.node(:,1),ScalpSurfaceMesh.node(:,2),ScalpSurfaceMesh.node(:,3),'.k')