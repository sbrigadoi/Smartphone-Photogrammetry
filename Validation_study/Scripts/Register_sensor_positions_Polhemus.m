clear all
close all
clc

addpath('../Utils')
load('../Other/phantom_positions.mat')

filename        = 'Operator4'; % Operator1/2/3/4
temp            = readtable(['../Polhemus/',filename,'.csv']);
pos             = table2array(temp(:,2:4));

idx_landmarks   = [153;154;288;134;135]; % Nasion, Inion, R, L, Cz
p1              = pos([1:5],:); % 
p2              = surface_points(idx_landmarks,:);

% Compute affine transformation between landmarks
T             = gen_xform_from_pts(p1,p2);
new_p1        = xform_apply(p1, T);

figure, plot3(p2(:,1), p2(:,2), p2(:,3), 'ob')
hold on, plot3(new_p1(:,1), new_p1(:,2), new_p1(:,3),'or')

% Apply affine transformation to all points
new_pos=xform_apply(pos, T);

% Find the correct label of each sensor
for i = 1:size(pos,1)
    dist                = sqrt(sum((new_pos(i,:)-surface_points).^2,2));
    [~, idx]            = min(dist);
    Correspondence_a(i)   = idx;
end

% Compute error 
error_a = sqrt(sum((new_pos-surface_points(Correspondence_a,:)).^2,2));

% Find the position corresponding to the second acquisition of the vertex
% (Cz)
idx     = find(Correspondence_a == 135);
idx(1)  = [];

% Recompute the affine transformation between landmarks using the second
% acquisition of the vertex position
p1              = pos([1:4 idx],:); % 
% Compute affine transformation between landmarks
T             = gen_xform_from_pts(p1,p2);
% Apply affine transformation to all points
new_pos=xform_apply(pos, T);

% Find the correct label of each sensor
for i = 1:size(pos,1)
    dist                = sqrt(sum((new_pos(i,:)-surface_points).^2,2));
    [~, idx]            = min(dist);
    Correspondence_b(i)   = idx;
end

% Compute error 
error_b = sqrt(sum((new_pos-surface_points(Correspondence_b,:)).^2,2));

if median(error_b) > median(error_a)
    error           = error_a;
    Correspondence  = Correspondence_a;
else
    error           = error_b;
    Correspondence  = Correspondence_b;
end
%%    

figure,
plot(error,'ob','MarkerFaceColor','b','MarkerSize',5)
xlim([0 290])
ylim([0 3.5])
xlabel('Points')
ylabel('Error [mm]')

save_file = ['../Results/',filename,'_error.mat'];
save(save_file,'error','Correspondence')




