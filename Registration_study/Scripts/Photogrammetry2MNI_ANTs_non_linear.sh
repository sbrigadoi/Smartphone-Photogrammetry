export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=96
export ANTSPATH= #Define ANTs path

subj= #name of the subject
T1_brain=../Data/MNI/MNI_no_nose.nii.gz
Atlas_brain= # nifti file with the subject's surface (obtained with the function Create_image_ANTs.m)
T1_brain_points=../Data/MNI/MNI_spheres.nii.gz
Atlas_brain_points= # nifti file with the landmarks positions (obtained with the function Create_image_ANTs.m)
sensor_position_file= # csv file with the sensor positions (obtained with the function Create_image_ANTs.m)
prefix_trasf=${subj}2MNI
path_result=../Data/${subj}/ANTs_non_linear/Sensor_positions/


cd ${path_result}
${ANTSPATH}/antsRegistration --verbose 1 --dimensionality 3 --float 1 \
--output ${prefix_trasf} \
-x [${T1_brain},${Atlas_brain}] \
-t rigid[0.1] \
-m PSE[${T1_brain_points},${Atlas_brain_points},1,1,1,600,5] \
--convergence [3000x1000x1000x1000,1e-8,10] \
--shrink-factors 8x4x2x1 \
--smoothing-sigmas 3x2x1x0vox \
-t SyN[0.25,6,0.0] \
-m PSE[${T1_brain},${Atlas_brain},1,1,1,600,5] \
--convergence [3000x1000x1000x1000,1e-8,10] \
--shrink-factors 8x4x2x1 \
--smoothing-sigmas 3x2x1x0vox \

${ANTSPATH}/antsApplyTransformsToPoints -d 3 -i ${sensor_position_file} -o ${prefix_trasf}.csv -t [${prefix_trasf}0GenericAffine.mat,1] -t ${prefix_trasf}1InverseWarp.nii.gz 
