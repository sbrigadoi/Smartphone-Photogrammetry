#!/bin/bash

# This script requires a previous brain extraction step
# To perform the brain extraction we use the following commands (T1.nii is the structural image)
# ${ANTSPATH}/N4BiasFieldCorrection -d 3 -i T1.nii -o N4T1.nii -s 2
# Matlab with nifti toolbox: reslice_nii(N4T1.nii,3DT1_RAS.nii)
# Mass software (https://www.nitrc.org/projects/cbica_mass/): mass -in 3DT1_RAS.nii -dest $folder_dataset -regs 15 -NOQ -MT 15 -regWt 0.3 -smooth 4 -v 
folder_dataset= # Define folder dataset
OUT_PREFIX_BRAIN=$folder_dataset/3DT1_RAS_2_MNI_ANTS_
OUT_WARPED_BRAIN=$folder_dataset/3DT1_RAS_T1_2_MNI_ANTS_Warped.nii.gz
OUT_INV_WARPED_BRAIN=$folder_dataset/3DT1_RAS_2_MNI_ANTS_InvWarped.nii.gz
template_brain=../Data/MNI/fmni_icbm152_t1_tal_nlin_asym_09c_brain.nii.gz
input_brain=$folder_dataset/3DT1_RAS_brain.nii.gz
export PATH=${PATH}: #Define path
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=96
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
ANTSPATH= #Define ANTs path
${ANTSPATH}/antsRegistration --verbose 1 --dimensionality 3 --float 0 --collapse-output-transforms 1 \
                    --output [$OUT_PREFIX_BRAIN,$OUT_WARPED_BRAIN,$OUT_INV_WARPED_BRAIN] \
                    --interpolation Linear --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] \
                    --initial-moving-transform [$template_brain,$input_brain,0] \
                    --transform Rigid[0.1] \
                                            --metric MI[$template_brain,$input_brain,1,32,Regular,0.25] \
                                            --convergence [500x500x500x500x500,1e-6,10] --shrink-factors 12x8x4x2x1 --smoothing-sigmas 4x3x2x1x0vox \
                    --transform Affine[0.1]  \
                                            --metric MI[$template_brain,$input_brain,1,32,Regular,0.25] \
                                            --convergence [500x500x500x500x500,1e-6,10] --shrink-factors 12x8x4x2x1 --smoothing-sigmas 4x3x2x1x0vox \
                    --transform SyN[0.1,3,0] \
                                            --metric CC[$template_brain,$input_brain,1,6] \
                                            --convergence [500x500x200x100x100,1e-6,10] --shrink-factors 12x8x4x2x1 --smoothing-sigmas 4x3x2x1x0vox
${ANTSPATH}/antsApplyTransforms --verbose 1 --dimensionality 3 --float 0 \
					--input $folder_dataset/sensor_positions2T1.nii  # Registration_study/Data/SubjX/Photogrammetry
					--reference-image $template \
					--output $folder_dataset/SensorPos2MNI.nii \  # Registration_study/Data/SubjX/T1
					--t $folder_dataset/3DT1_RAS_2_MNI_ANTS_1Warp.nii.gz \
					--t $folder_dataset/3DT1_RAS_2_MNI_ANTS_0GenericAffine.mat \
					-n NearestNeighbor \
