clear all
clc

ptCloud = pcread('../Meshes/Asus.ply');

cc      = ptCloud.Color;
loc     = double(ptCloud.Location);
col     = double(ptCloud.Color); 
hsv     = rgb2hsv(double(cc)./255);

%% 10-10 Positions

idx         = find(hsv(:,1)>0.2 & hsv(:,1)<0.5 & hsv(:,2)>0.07 & hsv(:,2)<0.8 & hsv(:,3)>0.3);
loc_10_10   = loc(idx,:);
figure, pcshow(ptCloud)
hold on,plot3(loc_10_10(:,1),loc_10_10(:,2),loc_10_10(:,3),'.c')

temp1 = nan(size(loc_10_10,1));
for i = 1:size(loc_10_10,1)
    current_pt = loc_10_10(i,:);
    temp1(i,:) = sqrt(sum((loc_10_10-current_pt).^2,2));
end
count = 1;
done  = [];
for i = 1:size(temp1,1)
    if length(find(i==done))==0
        current_dist    = temp1(i,:);
        nn              = find(current_dist<0.2);
        new             = [];
        for j = 1:length(nn)
            new_dist    = temp1(nn(j),:);
            new_nn      = find(new_dist<0.2);
            new         = [new new_nn];
        end
        new             = unique(new);
        nn              = unique([new nn]);
        cl(count).ngb   = nn;
        done            = [done nn];
        count           = count+1;
        plot3(loc_10_10(nn,1),loc_10_10(nn,2),loc_10_10(nn,3),'.k')
        %pause
        hold on
    end
end
hold off

pos_10_10 = nan(length(cl),3);
for i = 1:length(pos_10_10)
    pp              = cl(i).ngb;
    pp_pos          = loc_10_10(pp,:);
    pos_10_10(i,:)  = mean(pp_pos);
end
figure
pcshow(ptCloud)
hold on, plot3(pos_10_10(:,1), pos_10_10(:,2), pos_10_10(:,3), '*c')

clear temp1 cl nn cluster_size


%% 10-20 Positions

idx         = find(hsv(:,1)>0.6 & hsv(:,2)>0.5 & hsv(:,3)>0.3);
loc_10_20   = loc(idx,:);
figure, pcshow(ptCloud)
hold on, plot3(loc_10_20(:,1), loc_10_20(:,2), loc_10_20(:,3),'.c')

temp1 = nan(size(loc_10_20,1));
for i = 1:size(loc_10_20,1)
    current_pt = loc_10_20(i,:);
    temp1(i,:) = sqrt(sum((loc_10_20-current_pt).^2,2));
end
count = 1;
done  = [];
for i = 1:size(temp1,1)
    if length(find(i==done))==0
        current_dist    = temp1(i,:);
        nn              = find(current_dist<0.2);
        new             = [];
        for j = 1:length(nn)
            new_dist    = temp1(nn(j),:);
            new_nn      = find(new_dist<0.2);
            new         = [new new_nn];
        end
        new             = unique(new);
        nn              = unique([new nn]);
        cl(count).ngb   = nn;
        done            = [done nn];
        count           = count+1;
        plot3(loc_10_20(nn,1),loc_10_20(nn,2),loc_10_20(nn,3),'*k','MarkerSize',5)
        %pause
        hold on
    end
end
hold off

pos_10_20 = nan(length(cl),3);
for i = 1:length(pos_10_20)
    pp      = cl(i).ngb;
    pp_pos  = loc_10_20(pp,:);
    if size(pp_pos,1)>1
        pos_10_20(i,:)=mean(pp_pos);
    else
        pos_10_20(i,:)=pp_pos;
    end
end
figure
pcshow(ptCloud)
hold on, plot3(pos_10_20(:,1), pos_10_20(:,2), pos_10_20(:,3), '*c')

clear temp1 cl nn cluster_size

%% 10-5 Positions

idx         = find(hsv(:,1)>0.3 & hsv(:,1)<0.83);
loc_10_5     = loc(idx,:);
figure, pcshow(ptCloud)
hold on, plot3(loc_10_5(:,1), loc_10_5(:,2), loc_10_5(:,3), '.c')

temp1 = nan(size(loc_10_5,1));
for i = 1:size(loc_10_5,1)
    current_pt = loc_10_5(i,:);
    temp1(i,:) = sqrt(sum((loc_10_5-current_pt).^2,2));
end
count = 1;
done  = [];
for i = 1:size(temp1,1)
    if length(find(i==done))==0
        current_dist    = temp1(i,:);
        nn              = find(current_dist<0.15);
        new             = [];
        for j = 1:length(nn)
            new_dist    = temp1(nn(j),:);
            new_nn      = find(new_dist<0.15);
            new         = [new new_nn];
        end
        new             = unique(new);
        nn              = unique([new nn]);
        cl(count).ngb   = nn;
        done            = [done nn];
        count           = count+1;
        plot3(loc_10_5(nn,1), loc_10_5(nn,2), loc_10_5(nn,3),'.k')
        %pause
        hold on
    end
end
hold off

pos_10_5 = nan(length(cl),3);
for i  = 1:length(pos_10_5)
    pp              = cl(i).ngb;
    pp_pos          = loc_10_5(pp,:);
    if size(pp_pos,1)>1
        pos_10_5(i,:)=mean(pp_pos);
    else
        pos_10_5(i,:)=pp_pos;
    end
end

for i=1:size(pos_10_5,1)
    curr_dist   = sqrt(sum(((pos_10_5(i,:)-pos_10_10).^2),2));
    min_dist(i) = min(curr_dist);
end
figure, plot(min_dist,'o')
idx2=find(min_dist<0.1);
pos_10_5(idx2,:)=[];

figure
pcshow(ptCloud)
hold on, plot3(pos_10_5(:,1), pos_10_5(:,2), pos_10_5(:,3), '*c')

clear temp1 cl nn cluster_size

%% Landmark positions

idx     = find(hsv(:,1)>0.1 & hsv(:,1)<0.2 & hsv(:,2)>0.3);
loc_lm  = loc(idx,:);
figure, pcshow(ptCloud)
hold on, plot3(loc_lm(:,1), loc_lm(:,2), loc_lm(:,3), '.k')

temp1 = nan(size(loc_lm,1));
for i = 1:size(loc_lm,1)
    current_pt = loc_lm(i,:);
    temp1(i,:) = sqrt(sum((loc_lm-current_pt).^2,2));
end
count = 1;
done  = [];

for i = 1:size(temp1,1)
    if length(find(i==done))==0
        current_dist    = temp1(i,:);
        nn              = find(current_dist<0.3);
        new             = [];
        for j=1:length(nn)
            new_dist    = temp1(nn(j),:);
            new_nn      = find(new_dist<0.3);
            new         = [new new_nn];
        end
        new             = unique(new);
        nn              = unique([new nn]);
        cl(count).ngb   = nn;
        done            = [done nn];
        count           = count+1;
        plot3(loc_lm(nn,1),loc_lm(nn,2),loc_lm(nn,3),'.k')
%         pause
        hold on
    end
end
hold off

pos_lm = nan(length(cl),3);
for i = 1:length(pos_lm)
    pp          = cl(i).ngb;
    pp_pos      = loc_lm(pp,:);
    pos_lm(i,:) = mean(pp_pos);
end

figure
pcshow(ptCloud)
hold on, plot3(pos_lm(:,1),pos_lm(:,2),pos_lm(:,3),'*c')
%%
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

save_file = '../Results/Asus_pos.mat';
save(save_file, 'pos_lm', 'pos_10_10', 'pos_10_20', 'pos_10_5')



