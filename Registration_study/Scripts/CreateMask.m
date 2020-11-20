function [maskOriginal,maskCut] = createMask(pathname)
% This function creates the mask. The original MRI (T1 and T2)
% was segmented with SPM12 which outputs 6 tissue segmentation mask: c1 for
% the GM, c2 for the WM, c3 for the CSF, c4 for the skull, c5 for the skin
% and c6 for air. 
%
% Input:  pathname         pathname to the directory of the subject, which should
%                          contain the tissue probability maps from SPM12.
%
% Output: maskOriginal     segmented mask 
%         maskCut          segmented mask with 0 where neck and lower part
%                          of the head are present
%
% Written by S. Brigadoi 01/2017, updated 08/2018

% Load probability maps for the 6 tissue types and normalize to the maximum
% to get maps with values from 0 to 1.

cd pathname

pathnameSave = fullfile(pathname,'MaskMesh');

GM_prob  				= load_untouch_nii('c1N43DT1.nii');
GM_prob 				= GM_prob.img/max(GM_prob.img(:));
WM_prob  				= load_untouch_nii('c2N43DT1.nii');
WM_prob 				= WM_prob.img/max(WM_prob.img(:));
CSF_prob 				= load_untouch_nii('c3N43DT1.nii');
CSF_prob 				= CSF_prob.img/max(CSF_prob.img(:));
skull_prob 				= load_untouch_nii('c4N43DT1.nii');
temp1					= skull_prob.img;
temp1(find(temp1<245))	=	0;
skull_prob.img 			= temp1;
skull_prob 				= skull_prob.img/max(skull_prob.img(:));
scalp_prob 				= load_untouch_nii('c5N43DT1.nii');
temp1 					= scalp_prob.img;
temp1(find(temp1<252)) 	= 0;
scalp_prob.img 			= temp1;
scalp_prob 				= scalp_prob.img/max(scalp_prob.img(:));
air_prob 				= load_untouch_nii('c6N43DT1.nii');
temp1 					= air_prob.img;
temp1(find(temp1>20)) 	= 255;
air_prob.img 			= temp1;
air_prob 				= air_prob.img/max(air_prob.img(:));


% Plot to check number assigned to tissue by SPM is correct
figure;
subplot(321)
imagesc(rot90(GM_prob(:,:,round(end/2)+40)))
axis square
title('GM')
subplot(322)
imagesc(rot90(WM_prob(:,:,round(end/2)+40)))
axis square
title('WM')
subplot(323)
imagesc(rot90(CSF_prob(:,:,round(end/2)+40)))
axis square
title('CSF')
subplot(324)
imagesc(rot90(skull_prob(:,:,round(end/2)+40)))
axis square
title('Skull')
subplot(325)
imagesc(rot90(scalp_prob(:,:,round(end/2)+40)))
axis square
title('Soft tissue')
subplot(326)
imagesc(rot90(air_prob(:,:,round(end/2)+40)))
axis square
title('Air')
set(gcf,'PaperPositionMode','auto','Position',[560 85 1085 863])
%%
% Mask dimensions
dim = size(GM_prob);

% Compute segmentation mask and threshold it to get rid of external noise,
% thus creating a scalp mask. Air is not contained because the external
% space is segmented as air as well. 
mask_tot 					= WM_prob + CSF_prob + GM_prob + scalp_prob + skull_prob;
scalpMask 					= zeros(dim);
scalpMask(mask_tot > 0.99) 	= 1; 
for i = 1:size(scalpMask,3)
    scalpMask(:,:,i) = imfill(squeeze(scalpMask(:,:,i)),8,'holes');
%     figure,subplot(1,2,1),imagesc(squeeze(scalpMask1(:,i,:)))
%     subplot(1,2,2),imagesc(squeeze(scalpMask(:,i,:)))
%     pause
end
for i = 1:size(scalpMask,2)
    scalpMask(:,i,:) = imfill(squeeze(scalpMask(:,i,:)),8,'holes');
%     figure,subplot(1,2,1),imagesc(squeeze(scalpMask1(:,i,:)))
%     subplot(1,2,2),imagesc(squeeze(scalpMask(:,i,:)))
%     pause
end
for i = 1:size(scalpMask,1)
    scalpMask(i,:,:) = imfill(squeeze(scalpMask(i,:,:)),8,'holes');
%     figure,subplot(1,2,1),imagesc(squeeze(scalpMask1(:,i,:)))
%     subplot(1,2,2),imagesc(squeeze(scalpMask(:,i,:)))
%     pause
end
scalpMaskNew = scalpMask;

% Mask probability masks with scalp mask
GM_prob_brainMasked 					= zeros(dim);
GM_prob_brainMasked(scalpMaskNew==1) 	= GM_prob(scalpMaskNew == 1);

WM_prob_brainMasked 					= zeros(dim);
WM_prob_brainMasked(scalpMaskNew==1) 	= WM_prob(scalpMaskNew == 1);

CSF_prob_brainMasked 					= zeros(dim);
CSF_prob_brainMasked(scalpMaskNew==1) 	= CSF_prob(scalpMaskNew == 1);

scalp_prob_brainMasked 					= zeros(dim);
scalp_prob_brainMasked(scalpMaskNew==1) = scalp_prob(scalpMaskNew == 1);

skull_prob_brainMasked 					= zeros(dim);
skull_prob_brainMasked(scalpMaskNew==1) = skull_prob(scalpMaskNew == 1);

air_prob_brainMasked 					= zeros(dim);
air_prob_brainMasked(scalpMaskNew==1) 	= air_prob(scalpMaskNew == 1);

% Assign each voxel to the tissue with highest probability in that voxel
brain_tissues_t(:,:,:,1) = scalp_prob_brainMasked;
brain_tissues_t(:,:,:,2) = skull_prob_brainMasked;
brain_tissues_t(:,:,:,3) = CSF_prob_brainMasked;
brain_tissues_t(:,:,:,4) = GM_prob_brainMasked;
brain_tissues_t(:,:,:,5) = WM_prob_brainMasked;
brain_tissues_t(:,:,:,6) = air_prob_brainMasked;

[~,segm_prob] 	= max(brain_tissues_t,[],4);
segm_prob 		= segm_prob .*scalpMaskNew;
maskOriginal = segm_prob;


tiss_type = {'soft tissue','skull','CSF','GM','WM','air'};

% Plot to check segmentation mask
figure;
subplot(131)
imagesc(rot90(maskOriginal(:,:,round(end/2))));
axis equal tight
subplot(132)
imagesc(rot90(squeeze(maskOriginal(:,round(end/2),:))));
axis equal tight
subplot(133)
imagesc((rot90(squeeze(maskOriginal(round(end/2),:,:)))));
axis equal tight
set(gcf,'PaperPositionMode','auto','Position',[560 85 1085 863])

% Write the mask volume to the .mat file and the tissue type to
% the tissue type file
save('maskOriginal.mat','maskOriginal');
save('mask_tiss_type.txt', 'tiss_type');

% Remove voxels of neck and lower part of the head (set to 0)
maskCut = maskOriginal;
maskCut(:,1:80,:) = 0;

% Plot to check segmentation mask
figure;
subplot(131)
imagesc(rot90(maskCut(:,:,round(end/2))));
axis equal tight
subplot(132)
imagesc(rot90(squeeze(maskCut(:,round(end/2),:))));
axis equal tight
subplot(133)
imagesc((rot90(squeeze(maskCut(round(end/2),:,:)))));
axis equal tight
set(gcf,'PaperPositionMode','auto','Position',[560 85 1085 863])

% Write the mask volume to the .mat file 
save('maskCut.mat','maskCut');

% Save the mask to .nii file

T1 		= load_untouch_nii('N43DT1.nii');

T1.img 	= maskCut;
save_untouch_nii(T1,'temp.nii')
%%
reslice_nii('temp.nii','mask.nii',[],[],[],2)
a 	    = load_untouch_nii('mask.nii');
maskCut = double(a.img);
save maskCut maskCut

