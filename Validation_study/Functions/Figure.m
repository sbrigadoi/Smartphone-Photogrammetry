% This script creates Figure 5 of the manuscript

load('Validation_study/Other/phantom_positions.mat')

all_points      = 1:288;
files           = dir('Validation_study/Results/*error.mat');
error_matrix    = nan(length(all_points),length(files));

for i = 1:length(files)
    load(fullfile(files(i).folder,files(i).name));
    error_matrix(Correspondence,i) = error;
end

% With Polhemus we did not collect the lowest two rows of sensors. These
% positions were not even considered with the photogrammetry-based approaches.

idx_lowest_rows             = find((sum(isnan(error_matrix),2)));
error_matrix(idx_lowest_rows,:)    = [];  

for i = 1:length(files)
    labels{i} = files(i).name(1:end-10);
end

Photog          = {'Asus','Samsung','OnePlus','iPhone'};
Polhem          = {'Operator1','Operator2','Operator3','Operator4'};
for i = 1:length(Photog)
    idx_Pho(i) = find(strcmp(Photog{i},labels));
end
for i = 1:length(Polhem)
    idx_Pol(i) = find(strcmp(Polhem{i},labels));
end

x = error_matrix(:,[idx_Pho, idx_Pol]);

labels  = {'Asus','Samsung','OnePlus','iPhone','Op 1','Op 2','Op 3','Op 4'};
temp    =   round(linspace(0, 255,10));
col     = repmat(temp',2,3)/255;
group   = [ones(1,size(x,1)) 2*ones(1,size(x,1)) 3*ones(1,size(x,1)) ...
    4*ones(1,size(x,1))  5*ones(1,size(x,1))  6*ones(1,size(x,1)) ...
     7*ones(1,size(x,1))  8*ones(1,size(x,1))];
positions = [1 1.25 1.5 1.75 2.75 3 3.25 3.5];

figure,boxplot(x,group, 'positions', positions,'outliersize',4,'symbol','.k'); 

h = findobj(gca,'Tag','Median'); 
set(h,'Color','k'); 
xlim([0 4.5])
ylim([0 4])
ylabel('Error [mm]')
set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8))])
set(gca,'xticklabel',{'Photogrammetry','Polhemus'})
h = gca; h.XAxis.TickLength = [0 0];

h = findobj(gca,'Tag','Box');
for j = 1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),col(j),'FaceAlpha',0.5);
end
colormap(jet)
c = get(gca, 'Children');

hleg1 = legend(c(1:8),labels,'Units','Pixels','FontSize',8);
set(hleg1,'NumColumns',2,'box','off','Position',[190 320 1 1]);
txt = 'Photogrammetry';
text(35,310,txt,'FontSize',8,'Units','Pixels')
txt = 'Polhemus';
text(132,310,txt,'FontSize',8,'Units','Pixels')
