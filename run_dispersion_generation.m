% run generation

% N: 12, 23, 31, 78
clc, clear all , close all

N = [5, 10, 30, 50, 70]; % # of clusters
rd = 150*ones(1,5);
for i = 1:length(N)
    img_para = descriptor_recon_smooth(1000, 0.01, N(i), rd(i), 1);
    % plot
    figure()
%     plot(img_para(:,1),img_para(:,2),'.g')
    hold on
    h = ellipse(img_para(:,3), img_para(:,4), img_para(:,5)*180/3.1416, img_para(:,1),img_para(:,2),'r');
    axis equal
    save(['structure_output_N_',num2str(N(i))], 'img_para')
end

