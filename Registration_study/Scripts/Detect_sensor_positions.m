clear all
close all
clc

Subj    = '';
file    = ''; %Load the mesh of the subject derived from the photogrammetry-based method. Mesh should be in .ply format

ptCloud = pcread(file);
figure,pcshow(ptCloud)

cc      = ptCloud.Color;
loc     = double(ptCloud.Location); 
col     = double(ptCloud.Color); 
hsv     = rgb2hsv(double(cc)./255);

thresh  = 0.2;
thresh1 = 26;
thresh2 = 300;

%% Landmark positions identification
aa      = find(hsv(:,1)>0.5 & hsv(:,1)<0.575 & hsv(:,3)>0.4);
loc_lm  = loc(aa,:);
figure, pcshow(ptCloud)
hold on,plot3(loc_lm(:,1),loc_lm(:,2),loc_lm(:,3),'.c')
hold off

% figure
temp1 = nan(size(loc_lm,1));
for i = 1:size(loc_lm,1)
    current_pt = loc_lm(i,:);
    temp1(i,:) = sqrt(sum((loc_lm-current_pt).^2,2));
end
count   = 1;
done    = [];
for i = 1:size(temp1,1)
    if length(find(i==done))==0
        current_dist = temp1(i,:);
        nn           = find(current_dist<thresh);
        new          = [];
        for j = 1:length(nn)
            new_dist = temp1(nn(j),:);
            new_nn   = find(new_dist<thresh);
            new      = [new new_nn];
        end
        new             = unique(new);
        nn              = unique([new nn]);
        cl(count).ngb   = nn;
        done            = [done nn];
        count           = count+1;
%         plot3(loc_lm(nn,1),loc_lm(nn,2),loc_lm(nn,3),'.')
%         hold on
    end
end
% hold off

for i = 1:length(cl)
    ll(i)  = length(cl(i).ngb);
end
cl(find(ll<thresh1 | ll>thresh2)) = [];

figure
pcshow(ptCloud)
hold on
for i = 1:length(cl)
    pp          = cl(i).ngb;
    pp_pos      = loc_lm(pp,:);
    pos_lm(i,:) = mean(pp_pos);
    plot3(pos_lm(i,1),pos_lm(i,2),pos_lm(i,3),'*c','MarkerSize',20)
    hold on
end
hold off 
% The order of the landmarks must be manually checked
% In this case the order is: right preauricular point,left preauricular point, inion, nasion

%% 10_20 positions identification
aa          = find(hsv(:,1)<0.04 & hsv(:,3)>0.78); 
loc_10_20   = loc(aa,:);
pcshow(ptCloud)
hold on,plot3(loc_10_20(:,1),loc_10_20(:,2),loc_10_20(:,3),'.c')

figure
temp1 = nan(size(loc_10_20,1));
for i = 1:size(loc_10_20,1)
    current_pt = loc_10_20(i,:);
    temp1(i,:) = sqrt(sum((loc_10_20-current_pt).^2,2));
end
count   = 1;
done    = [];
for i = 1:size(temp1,1)
    if length(find(i==done))==0
        current_dist    = temp1(i,:);
        nn              = find(current_dist<thresh);
        new             = [];
        for j = 1:length(nn)
            new_dist    = temp1(nn(j),:);
            new_nn      = find(new_dist<thresh);
            new         = [new new_nn];
        end
        new             = unique(new);
        nn              = unique([new nn]);
        cl1(count).ngb  = nn;
        done            = [done nn];
        count           = count+1;
        plot3(loc_10_20(nn,1),loc_10_20(nn,2),loc_10_20(nn,3),'.c')
        hold on
    end
end
hold off

for i = 1:length(cl1);
    ll1(i)=length(cl1(i).ngb);
end
cl1(find(ll1<thresh1 | ll1>thresh2)) = [];

figure
pcshow(ptCloud)
hold on
for i = 1:length(cl1)
    pp              = cl1(i).ngb;
    pp_pos          = loc_10_20(pp,:);
    pos_10_20(i,:)  = mean(pp_pos);
    plot3(pos_10_20(i,1),pos_10_20(i,2),pos_10_20(i,3),'*c','MarkerSize',20)        
    hold on
end
% Cz is at position 19

%% 10-10 positions identification
aa          = find(hsv(:,1)>0.1 & hsv(:,1)<0.5 & hsv(:,2)>0.3 & hsv(:,2)<0.7 & hsv(:,3)>0.4);
loc_10_10   = loc(aa,:);
pcshow(ptCloud)
hold on,plot3(loc_10_10(:,1),loc_10_10(:,2),loc_10_10(:,3),'.c')

figure
% pcshow(ptCloud)
% hold on
temp1 = nan(size(loc_10_10,1));
for i = 1:size(loc_10_10,1)
    current_pt = loc_10_10(i,:);
    temp1(i,:) = sqrt(sum((loc_10_10-current_pt).^2,2));
end
count   = 1;
done    = [];
for i=1:size(temp1,1)
    if length(find(i==done))==0
        current_dist    = temp1(i,:);
        nn              = find(current_dist<thresh);
        new             = [];
        for j=1:length(nn)
            new_dist    = temp1(nn(j),:);
            new_nn      = find(new_dist<thresh);
            new         = [new new_nn];
        end
        new             = unique(new);
        nn              = unique([new nn]);
        cl2(count).ngb  = nn;
        done            = [done nn];
        count           = count+1;
        plot3(loc_10_10(nn,1),loc_10_10(nn,2),loc_10_10(nn,3),'.')
        hold on
    end
end
hold off

for i=1:length(cl2);
    ll2(i) = length(cl2(i).ngb);
end
cl2(find(ll2<thresh1 | ll2>thresh2)) = [];

figure
pcshow(ptCloud)
hold on
for i=1:length(cl2)
    pp              = cl2(i).ngb;
    pp_pos          = loc_10_10(pp,:);
    pos_10_10(i,:)  = mean(pp_pos);
    plot3(pos_10_10(i,1),pos_10_10(i,2),pos_10_10(i,3),'*c','MarkerSize',20)
    hold on
end

%% 
pos_cz          = pos_10_20(19,:); 
pos_10_20(19,:) = []; % Delete Cz position from the set of positions of the 10_20 system
% oopy Cz in lm(5).pos
lm_pos(1).name  = 'right';
lm_pos(1).pos   = pos_lm(1,:);
lm_pos(2).name  = 'left';
lm_pos(2).pos   = pos_lm(2,:);
lm_pos(3).name  = 'inion';
lm_pos(3).pos   = pos_lm(3,:);
lm_pos(4).name  = 'nasion';
lm_pos(4).pos   = pos_lm(4,:);
lm_pos(5).name  = 'cz';
lm_pos(5).pos   = pos_cz;

%%  
load(fullfile('..','Other','phantom_positions.mat'))
[num,txt,raw]   = xlsread(fullfile('..','Other','name_positions.xlsx'));
label_ref       = {txt{2:289,1}};
label_landmarks = {'right','left','inion','nasion','cz'};
label_10_20     = {'Fz','Pz','Fp2','F4','C4','P4','O2','F8','T8','P8','Fp1','F3','C3','P3','O1','F7','T7','P7'};
label_10_10     = {'Fpz','AF8','FT8','TP8','PO8','Oz','PO7','TP7','FT7','AF7', ...
    'AFz','AF4','F6','FC6','C6','CP6','P6','PO4','POz','PO3','P5','CP5','C5','FC5','F5','AF3', ...
    'F2','FC4','CP4','P2','P1','CP2','CP3','FC3','F1','CPz','CP1', ...
    'FCz','FC2','C2','C1','FC1'}; %we deleted the positions in the last row because these positions are not in the phantom 'F10','FT10','TP10','P10','PO10','O10','O9','PO9','P9','TP9','FT9','F9'

for i = 1:length(label_landmarks)
    idx_landmarks(i) = find(strcmp(label_landmarks{i},label_ref));
end
for i = 1:length(label_10_20)
    idx_10_20(i)    = find(strcmp(label_10_20{i},label_ref));
end
for i = 1:length(label_10_10)
    idx_10_10(i)    = find(strcmp(label_10_10{i},label_ref));
end

%% The following steps aim to assign the correct label to each detected position

addpath(genpath(fullfile('..','Utils')))

for i = 1:length(idx_landmarks)
    p1(i,:) = lm_pos(i).pos;
end

p2      = surface_points(idx_landmarks,:);
T       = gen_xform_from_pts(p1,p2);
% new_p1  = xform_apply(p1, T);
% figure,plot3(p2(:,1),p2(:,2),p2(:,3),'ob')
% hold on,plot3(new_p1(:,1),new_p1(:,2),new_p1(:,3),'or')

temp    = xform_apply(pos_10_10, T);
temp2   = surface_points(idx_10_10,:);
clear opt Transform
opt.method      = 'rigid'; 
opt.viz         = 1; 
opt.normalize   = 1; 
opt.scale       = 1;
opt.rot         = 1;         
opt.corresp     = 0;
opt.tol         = 1e-8;
opt.corresp     = 1;
Y               = temp; 
X               = temp2;
[Transform, C]  = cpd_register(X,Y,opt);
ee              = sum((Transform.Y-X(C,:)).^2,2);figure,plot(ee,'*')
last_row        = find(ee > 500);
figure, plot3(Transform.Y(:,1),Transform.Y(:,2),Transform.Y(:,3),'*b'), hold on,
plot3(Transform.Y(last_row,1),Transform.Y(last_row,2),Transform.Y(last_row,3),'or')
legend('Sensor positions','Sensor belonging to the last row')

pos_10_10(last_row,:)  = [];
Y(last_row,:)          = [];
[Transform, C]         = cpd_register(X,Y,opt);
l10                    = {label_ref{idx_10_10}};
l10                    = {l10{C}};

for i = 1:length(l10)
    s10_10_pos(i).name = l10{i};
    s10_10_pos(i).pos  = pos_10_10(i,:);

end
    
%%
temp            = xform_apply(pos_10_20, T);
temp2           = surface_points(idx_10_20,:);
clear opt Transform
opt.method      = 'rigid'; 
opt.viz         = 1;
opt.normalize   = 1;
opt.scale       = 1;
opt.rot         = 1;         
opt.corresp     = 0;
opt.tol         = 1e-8; 
opt.corresp     = 1;
Y               = temp; 
X               = temp2;
[Transform, C]  = cpd_register(X,Y,opt);
l20             = {label_ref{idx_10_20}};
l20             = {l20{C}};

for i = 1:length(l20)
    s10_20_pos(i).name = l20{i};
    s10_20_pos(i).pos  = pos_10_20(i,:);
end

%%
save_file = fullfile('..','Data',Subj,'Photogrammetry','sensor_positions.mat');
save(save_file,'lm_pos','s10_10_pos','s10_20_pos')

