close all, clc

N = [12, 23, 31, 78]; % # of clusters
for i = 1:length(N)
load (['structure_output_N_',num2str(N(i))])
    figure()

    h = ellipse(img_para(:,3), img_para(:,4), img_para(:,5)*180/3.1416, img_para(:,1),img_para(:,2),'k');
    axis equal
end