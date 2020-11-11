clear all
close all
clc

addpath('../Utils')
load('../Other/phantom_positions.mat')
load('../Results/OnePlus_pos.mat')
ptCloud         = pcread('../Meshes/OnePlus.ply');
[num,txt,raw]   = xlsread('../Other/name_positions.xlsx');

landmarks = {'L3','L4','Inion','Nasion','Cz'};

label_ref = {txt{2:289,1}};
for i = 1:length(landmarks)
    idx_landmarks(i) = find(strcmp(landmarks{i}, label_ref));
end

% Order_pos_lm: Nasion,R,S,In,Cz
p1(1,:) = pos_lm(2,:);
p1(2,:) = pos_lm(3,:);
p1(3,:) = pos_lm(4,:);
p1(4,:) = pos_lm(1,:);
p1(5,:) = pos_lm(5,:);
p2      = surface_points(idx_landmarks,:);

%Compute affine transformation between landmarks
T       = gen_xform_from_pts(p1,p2);
new_p1  = xform_apply(p1, T);

%%
figure, plot3(p2(:,1), p2(:,2), p2(:,3), 'ob')
hold on, plot3(new_p1(:,1), new_p1(:,2), new_p1(:,3), 'or')


% Apply affine transformation to all points
pos     = [pos_10_10; pos_10_20; pos_10_5; p1];
new_pos = xform_apply(pos, T);

% Find the correct label of each sensor
for i = 1:size(pos,1)
    dist                = sqrt(sum((new_pos(i,:)-surface_points).^2,2));
    [~, idx]            = min(dist);
     Correspondence(i)  = idx;
end

% Compute error 
error = sqrt(sum((new_pos-surface_points(Correspondence,:)).^2,2));

figure,
plot(error,'ob','MarkerFaceColor','b','MarkerSize',5)
xlim([0 290])
ylim([0 3.5])
xlabel('Points')
ylabel('Error [mm]')

save_file = '../Results/OnePlus_error.mat';
save(save_file,'error','Correspondence')

