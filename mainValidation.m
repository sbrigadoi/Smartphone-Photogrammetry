% This script performs the validation study 

% Mobile 1: Asus
% Mobile 2: Samsung
% Mobile 3: OnePlus
% Mobile 4: iPhone

addpath(fullfile('Validation_study','Functions'))

mobileName = {'Asus','Samsung','OnePlus','iPhone'};
distThresh = [0.15 0.3 0.2 0.3];
filenamePol = {'Operator1','Operator2','Operator3','Operator4'};

% Load ground truth (phantom positions and their labels)
load(fullfile('Validation_study','Other','phantom_positions.mat'))
[~,txt,~] = xlsread(fullfile('Validation_study','Other','name_positions.xlsx'));
landmarks = {'L3','L4','Inion','Nasion','Cz'};
label_ref = {txt{2:289,1}};
idx_landmarks = zeros(length(landmarks),1);
for i = 1:length(landmarks)
    idx_landmarks(i) = find(strcmp(landmarks{i}, label_ref));
end
p2 = surface_points(idx_landmarks,:);
     
% Photogrammetry
for iM = 1:length(mobileName)
    
    % Read mesh file created from the video recorded by the mobile
    ptCloud = pcread(fullfile('Validation_study','Meshes',[mobileName{iM} '.ply']);

    cc      = ptCloud.Color;
    loc     = double(ptCloud.Location);
    col     = double(ptCloud.Color);
    hsv     = rgb2hsv(double(cc)./255);
    
    % Detect 10-10 Positions (thresholds on colors depend on video quality)
    switch mobileName{iM}
        case 'Samsung'
            idx = find(hsv(:,1)>0.2 & hsv(:,1)<0.5 & hsv(:,2)>0.05 & hsv(:,2)<0.9);
        case 'OnePlus'
            idx = find(hsv(:,1)>0.2 & hsv(:,1)<0.5 & hsv(:,2)>0.11 & hsv(:,2)<0.9);
        case 'iPhone'
            idx = find(hsv(:,1)>0.3 & hsv(:,1)<0.5 & hsv(:,2)>0.3 & hsv(:,2)<0.7 & hsv(:,3)>0.3);
        case 'Asus'
            idx = find(hsv(:,1)>0.2 & hsv(:,1)<0.5 & hsv(:,2)>0.07 & hsv(:,2)<0.8 & hsv(:,3)>0.3);
    end
    loc_10_10 = loc(idx,:);
    
%     % Display detected 10-10 positions
%     figure;
%     pcshow(ptCloud)
%     hold on; 
%     plot3(loc_10_10(:,1), loc_10_10(:,2), loc_10_10(:,3), '.c')
    
    % Identify the real coordinate of each position
    pos_10_10 = identifyPoint(loc_10_10,distThresh(iM));
    
%     figure
%     pcshow(ptCloud)
%     hold on, plot3(pos_10_10(:,1), pos_10_10(:,2), pos_10_10(:,3), '*c')
        
    % Detect 10-20 Positions (thresholds on colors depend on video quality)
    switch mobileName{iM}
        case 'Samsung'
            idx = find(hsv(:,1)<0.05 & hsv(:,2)>0.37 & hsv(:,3)>0.3);
        case 'OnePlus'
            idx = find(hsv(:,1)>0.6 & hsv(:,2)>0.6 & hsv(:,3)>0.3);
        case 'iPhone'
            idx = find(hsv(:,1)<0.11 & hsv(:,2)>0.6  & hsv(:,3)>0.45);
        case 'Asus'
            idx = find(hsv(:,1)>0.6 & hsv(:,2)>0.5 & hsv(:,3)>0.3);
    end
    loc_10_20   = loc(idx,:);
    
%     % Display detected 10-20 positions
%     figure;
%     pcshow(ptCloud)
%     hold on;
%     plot3(loc_10_20(:,1), loc_10_20(:,2), loc_10_20(:,3), '.c')

    % Identify the real coordinate of each position
    pos_10_20 = identifyPoint(loc_10_20,distThresh(iM));
    
%     figure
%     pcshow(ptCloud)
%     hold on, plot3(pos_10_20(:,1), pos_10_20(:,2), pos_10_20(:,3),'*c')
    
    % Detect 10-5 Positions (thresholds on colors depend on video quality)
    switch mobileName{iM}
        case 'Samsung'
            idx = find(hsv(:,1)>0.45 & hsv(:,1)<0.9);
        case 'OnePlus'
            idx = find(hsv(:,1)>0.45 & hsv(:,1)<0.9 & hsv(:,2)>0.12);
        case 'iPhone'
            idx = find(hsv(:,1)>0.5 & hsv(:,1)<0.7 & hsv(:,2)>0.2);
        case 'Asus'
            idx = find(hsv(:,1)>0.3 & hsv(:,1)<0.83);
    end
    loc_10_5    = loc(idx,:);
    
%     % Display detected 10-5 positions
%     figure;
%     pcshow(ptCloud)
%     hold on;
%     plot3(loc_10_5(:,1), loc_10_5(:,2), loc_10_5(:,3), '.c')
    
    % Identify the real coordinate of each position
    pos_10_5 = identifyPoint(loc_10_5,distThresh(iM));
    
    if strcmpi(mobileName{iM},'Asus') % Further check for Asus
        for i=1:size(pos_10_5,1)
            curr_dist   = sqrt(sum(((pos_10_5(i,:)-pos_10_10).^2),2));
            min_dist(i) = min(curr_dist);
        end
        idx2=find(min_dist<0.1);
        pos_10_5(idx2,:)=[];
    end
    
%     figure
%     pcshow(ptCloud)
%     hold on, plot3(pos_10_5(:,1), pos_10_5(:,2), pos_10_5(:,3), '*c')
    
    % Detect landmark positions (thresholds on colors depend on video quality)
    switch mobileName{iM}
        case 'Samsung'
            idx = find(hsv(:,1)>0.05 & hsv(:,1)<0.3 & hsv(:,2)>0.45);
        case 'OnePlus'
            idx = find(hsv(:,1)>0.1 & hsv(:,1)<0.3 & hsv(:,2)>0.25 & hsv(:,3)>0.4);
        case 'iPhone'
            idx = find(hsv(:,1)>0.12 & hsv(:,1)<0.3 & hsv(:,2)>0.5);
        case 'Asus'
            idx = find(hsv(:,1)>0.1 & hsv(:,1)<0.2 & hsv(:,2)>0.3);
    end
    loc_lm  = loc(idx,:);
    
%     % Display detected landmark positions
%     figure;
%     pcshow(ptCloud)
%     hold on;
%     plot3(loc_lm(:,1), loc_lm(:,2), loc_lm(:,3), '.k')

    % Identify the real coordinate of each position
    pos_lm = identifyPoint(loc_lm,distThresh(iM));
     
%      figure
%      pcshow(ptCloud)
%      hold on, plot3(pos_lm(:,1),pos_lm(:,2),pos_lm(:,3),'*c')
     
     % Summary figure
     figure
     pcshow(ptCloud)
     hold on, plot3(pos_lm(:,1),pos_lm(:,2),pos_lm(:,3),'*y','MarkerSize',15,'LineWidth',2)
     hold on, plot3(pos_10_10(:,1),pos_10_10(:,2),pos_10_10(:,3),'*g','MarkerSize',15,'LineWidth',2)
     hold on, plot3(pos_10_20(:,1),pos_10_20(:,2),pos_10_20(:,3),'*r','MarkerSize',15,'LineWidth',2)
     hold on, plot3(pos_10_5(:,1),pos_10_5(:,2),pos_10_5(:,3),'*b','MarkerSize',15,'LineWidth',2)
     grid off
     set(gca,'xtick',[],'ytick',[],'ztick',[])
     axis off
     view([20,150])
     
     % Saving 
%      save_file = fullfile('Validation_study','Results',[mobileName{iM} '_pos.mat']);
%      save(save_file, 'pos_lm', 'pos_10_10', 'pos_10_20', 'pos_10_5')

     % REGISTRATION
     % Register the positions to the head template
     % Order_pos_lm: Inion,Cz,R,L,Nasion
     switch mobileName{iM}
        case 'Samsung'
            p1 = [pos_lm(3,:); pos_lm(4,:); pos_lm(1,:); pos_lm(5,:); pos_lm(2,:)];
        case 'OnePlus'
            p1 = [pos_lm(2,:); pos_lm(3,:); pos_lm(4,:); pos_lm(1,:); pos_lm(5,:)];
        case 'iPhone'
            p1 = [pos_lm(3,:); pos_lm(4,:); pos_lm(1,:); pos_lm(5,:); pos_lm(2,:)];
        case 'Asus'
            p1 = [pos_lm(3,:); pos_lm(2,:); pos_lm(5,:); pos_lm(1,:); pos_lm(4,:)];
     end

     % Compute affine transformation between landmarks
     T       = affineTransf(p1,p2);
     new_p1  = applyTransf(p1, T);

%      figure, plot3(p2(:,1), p2(:,2), p2(:,3), 'ob')
%      hold on, plot3(new_p1(:,1), new_p1(:,2), new_p1(:,3), 'or')

     % Apply affine transformation to all points
     pos     = [pos_10_10; pos_10_20; pos_10_5; p1];
     new_pos = applyTransf(pos, T);

     % Find the correct label of each sensor
     Correspondence = zeros(size(pos,1),1);
     for i = 1:size(pos,1)
         dist                = sqrt(sum((new_pos(i,:)-surface_points).^2,2));
         [~, idx]            = min(dist);
         Correspondence(i,1)  = idx;
     end
     
     % Compute error
     error = sqrt(sum((new_pos-surface_points(Correspondence,:)).^2,2));
     
     figure,
     plot(error,'ob','MarkerFaceColor','b','MarkerSize',5)
     xlim([0 290])
     ylim([0 3.5])
     xlabel('Points')
     ylabel('Error [mm]')

     % Saving
%      save_file = fullfile('Validation_study','Results',[mobileName{iM} '_error.mat']);
%      save(save_file,'error','Correspondence')
end

% Polhemus 
for iO = 1:length(filenamePol)
    
    % Read Polhemus file 
    temp            = readtable(fullfile('Validation_study','Polhemus',[filenamePol{iO},'.csv']));
    pos             = table2array(temp(:,2:4));
    
    idx_landmarks   = [153;154;288;134;135]; % Nasion, Inion, R, L, Cz
    p1              = pos(1:5,:); 
    p2              = surface_points(idx_landmarks,:);

    % Compute affine transformation between landmarks
    T             = affineTransf(p1,p2);
    new_p1        = applyTransf(p1,T);

    % figure, plot3(p2(:,1), p2(:,2), p2(:,3), 'ob')
    % hold on, plot3(new_p1(:,1), new_p1(:,2), new_p1(:,3),'or')

    % Apply affine transformation to all points
    new_pos=applyTransf(pos,T);

    % Find the correct label of each sensor
    Correspondence_a = zeros(size(pos,1),1);
    for i = 1:size(pos,1)
        dist                = sqrt(sum((new_pos(i,:)-surface_points).^2,2));
        [~, idx]            = min(dist);
        Correspondence_a(i,1)   = idx;
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
    T             = affineTransf(p1,p2);
    % Apply affine transformation to all points
    new_pos = applyTransf(pos, T);

    % Find the correct label of each sensor
    Correspondence_b = zeros(size(pos,1),1);
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
    
    figure,
    plot(error,'ob','MarkerFaceColor','b','MarkerSize',5)
    xlim([0 290])
    ylim([0 3.5])
    xlabel('Points')
    ylabel('Error [mm]')
    
    % Saving
%     save_file = ['Validation_study/Results/',filenamePol{iO},'_errorProva.mat'];
%     save(save_file,'error','Correspondence')
end




