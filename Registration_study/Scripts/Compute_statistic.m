clear all
close all
clc

addpath(genpath(fullfile('..','Utils')))
cd (fullfile('..','Results'))
c_type = {'Affine_landmarks','CPD_affine','CPD_non_linear','ANTs_affine','ANTs_non_linear'};

for i = 1:length(c_type)
    load(['Sensor_pos_',c_type{i},'.mat'])
    values(1,:,i) = euclidean_dist(:);
    values(2,:,i) = geodesic_dist(:);
    load(['Fluence_act_pos_',c_type{i},'.mat'])
    values(3,:,i) = euclidean_dist(:);
    values(4,:,i) = geodesic_dist(:);
    load(['Fluence_distrib_dice_',c_type{i},'.mat'])
    values(5,:,i) = dice(:);
end

combos = combntns([1:5],2);
m_type = {'Eu_sens','Ge_sens','Eu_fl','Ge_fl','Dice_fl'};
for mm = 1:length(m_type)
    temp = squeeze(values(mm,:,:));
    for cc = 1:size(combos,1)
        [p,h,stat]       = signrank(temp(:,combos(cc,1)),temp(:,combos(cc,2)),'tail','both');
        p_mat(mm,cc)     = p;
        sign_test(mm,cc) = stat.zval;
    end
    fdr_mat(mm,:) = fdr_bh(p_mat(mm,:),0.05);
    T  = table({c_type{combos(:,1)}}',{c_type{combos(:,2)}}',fdr_mat(mm,:)',sign_test(mm,:)',p_mat(mm,:)');
    T.Properties.VariableNames{'Var3'} = 'H';
    T.Properties.VariableNames{'Var4'} = 'Sign';
    T.Properties.VariableNames{'Var5'} = 'p';
    eval([m_type{mm},'_table= T; '])
end





