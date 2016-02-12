function [ img_para ] = descriptor_recon_smooth( L, VF, N, ndt, Compact )
% -------------------------------------------------------------------------
% Descriptor-based reconstruction of microstructure
% the output is a structural parameter matrix for constructing
% smooth-profile images
%
% With minor revision, it can also generate pixelated images.
% 
% Inputs:
% L: side length of reconstructed image
% VF: volume fraction
% N: number of cluster
% ndt: nearest center distance of clusters
% Compact: aspect ratio
%
% Outputs:
% img_para: a matrix record the information of all clusters in the image.
% Each row is one cluster;
% Each column represent one structural parameters:
% [ Center X;  Center Y;  Major radius;  Minor radius;  Oriantation angle ]
% 
% FUNCTIONS USED:
% remove_single.m
% (find_edge_0.m, for pixelated)
% (find_edge_1.m, for pixelated)
%
% By Hongyi Xu, hongyixu2014@u.northwestern.edu
% -------------------------------------------------------------------------
img_para = zeros(N,5);

ps = 2;  % The size of patch (single cluster box). The side-length of the box = L/ps. Decrease ps when domain decrease or cluster size increase to avoid error
%% Generate a random microstructure
all = 1:1:(L*L);
all = all';
[temp1, temp2] = meshgrid(1:L,1:L);
temp1 = temp1(:);
temp2 = temp2(:);
all = [all, temp1, temp2];
clear temp1 temp2

xl = [];  yl = [];  % the list for x and y coordinates
img = zeros(L,L);
a = randperm(L*L); % 2500 = 50*50, the total number of the pixels.
for ii=1:1:N % N points are 1
    jj = a(ii);
    xl = [ xl; all(jj,2) ];
    yl = [ yl; all(jj,3) ];
    img( all(jj,2) , all(jj,3) )=1;  % img should still be kept. In random moving, it is convenient to use img to check overlap
end
temp_img = zeros(L+2,L+2);
temp_img(2:L+1,2:L+1) = img;  % add a 1 pixel wide "buff region"
img = temp_img;
xl = xl + 1;
yl = yl + 1;
clear ii jj a 
%% Adjust the nearest distances using Simulated Annealing
% Evaluate the nearest distances for the first time
cl = [xl,yl];  % Coordinate list
nd = [];
for ii = 1:1:N
    expand_c = repmat( cl(ii,:) , [N , 1] );
    distances = (cl - expand_c).^2;
    distances = sum( distances , 2);
    distances = distances.^0.5;
    distances = sort(distances);
    distances = distances( 2:length(distances) );
    nd = [nd; min(distances)];
end
% Start Simulated Annealing
disp('Now start to reconstruct the dispersion of cluster centers ...')
T = 100;
Pr = 0.5;
while T > 0  % T: change from 200 -> 0
    
    move_order = randperm(N);
    for ii = 1:1:N
        
        nd_old = nd;
        num = move_order(ii);
        x = cl(num,1);
        y = cl(num,2);
        while 1
            
            while 1
                x_move = 1 + ceil( -5 + 9*rand() );
                y_move = 1 + ceil( -5 + 9*rand() );
                if x+x_move >= 2 && x+x_move <= L+1  && y+y_move >= 2 && y+y_move <= L+1
                    break
                end
            end
%             if img( x+x_move-1 , y+y_move ) + img( x+x_move+1 , y+y_move ) + img( x+x_move , y+y_move-1 ) + img( x+x_move , y+y_move+1 ) + img( x+x_move , y+y_move ) == 0  % No adjacsent 1 pixels
            if sum( sum( img(x+x_move-1:x+x_move+1 , y+y_move-1:y+y_move+1) ) ) == 0  % or use the observed minimum nearest center distance  
                % Re-evaluate
                
                cl(num,:) = [ x+x_move , y+y_move ];
                nd = [];
                for jj = 1:1:N
                    expand_c = repmat( cl(jj,:) , [N , 1] );
                    distances = (cl - expand_c).^2;
                    distances = sum( distances , 2);
                    distances = distances.^0.5;
                    distances = sort(distances);
                    distances = distances( 2:length(distances) );
                    nd = [nd; min(distances)];
                end
                
%                 distances = (cl - expand_c).^2;
%                 distances = sum( distances , 2);
%                 distances = distances.^0.5;
%                 distances = sort(distances);
%                 distances = distances( 2:length(distances) );
%                 nd(num) = min(distances);

                % If better, accept; if worse, accept with a given Pr.
                if abs( mean(nd) - ndt ) < abs( mean(nd_old) - ndt )  % if better
                    img(x,y) = 0;
                    img(x+x_move, y+y_move) = 1;
%                     cl(num,1) = x+x_move;
%                     cl(num,2) = y+y_move;
                else  % if worse
                    flag = rand();
                    if flag < Pr  % accept the bad with certain Pr
                        img(x,y) = 0;
                        img(x+x_move, y+y_move) = 1;
%                         cl(num,1) = x+x_move;
%                         cl(num,2) = y+y_move;
                    else  % not accept
                        nd = nd_old;
                        cl(num,1) = x;
                        cl(num,2) = y;
                    end
                end
                break;
            end
            
        end
        
    end
%     [x y] = find(img == 1);
%     cll = [x,y];
    T = T - 1;
    Pr = Pr - (250-T)*Pr/250;
    disp(T);
    disp( mean(nd) );
    if abs( mean(nd) - ndt ) < 0.01
        break;
    end
    
end
img = img( 2:L+1 , 2:L+1 );
[x y] = find(img == 1);
cl = [x,y];
% Now we finish the nearest center distances. The only problem is that some
% center pixels may be too close to each other, even though they are not
% touched.
img_para(:,1:2) = cl;

%% Assign area and compactness to each cluster
disp('Now assigning the area and compactness values for each cluster...')
% % --------------- Gaussian distribution of cluster areas ------------------
% Area_m = L^2*VF/N;
% Area_v = Area_m/10;
% areas = randn([N,1]);
% areas = areas*Area_v + Area_m;
% % -------------------------------------------------------------------------

% % ---------------- Uniform distribution of cluster areas ------------------
% Area_m = L^2*VF/N;
% areas = rand([N,1]);
% areas = areas * Area_m/mean(areas);
% % -------------------------------------------------------------------------

% --------------- Exponential distribution of cluster areas ---------------
mu = L^2*VF/N;
areas = exprnd(mu,[N,1]);
% -------------------------------------------------------------------------

ng = find(areas < 0);  % ng: negative values
if ng>0
    for ii = 1:1:length(ng)
        areas( ng(ii) ) = 1;  % use 1 to replace the negative value. The min bias
    end
end

areas =  areas / ( sum( areas(:) )/L^2 / VF );

Comp_m = Compact;
Comp_v = Comp_m/10;
comps = randn([N,1]);
comps = comps*Comp_v + Comp_m;
ng = find(comps < 0);  % ng: negative values
if ng>0
    for ii = 1:1:length(ng)
        comps( ng(ii) ) = 1;  % use 1 to replace the negative value. The min bias
    end
end
%% generate single cluster images
SIS = {};  % all the generated single images
area_sum = 0;
for ii = 1:1:N
    si = zeros(L/ps,L/ps);  % si: single image
    A = areas(ii);
    elong = comps(ii);
    a = sqrt( A/pi/elong );
    b = a*elong;
    for x = (L/ps/2 - ceil(a) ):1:(L/ps/2 + ceil(a) )
        for y = (L/ps/2 - ceil(b) ):1:(L/ps/2 + ceil(b) )
            if ( (x - (L/ps)/2)/a )^2 + ( (y - (L/ps)/2)/b )^2 <= 1
                si(x,y) = 1;
            end
        end
    end
    
    % find the [x y] coordinates for the '1's in the binary image
    [x y] = find(si == 1);
    % adjust the (0,0)point to the center
    tc = [x y];  % temp coordinate
    tc = tc - (L/ps/2)-0.5;
    theta = 360*rand();  % Can change the tilt angle range. Full rotation range = 180
    
    img_para(ii,3) = a;
    img_para(ii,4) = b;
    img_para(ii,5) = theta;
    
    theta = theta/180*pi;
    tc = tc * [cos(theta)  -sin(theta); sin(theta)   cos(theta)];
    tc = tc + (L/ps/2)+0.5;
    tc = floor(tc);
    si = zeros(L/ps,L/ps);
    for jj = 1:1:size(tc,1)
        si( tc(jj,1) , tc(jj,2) ) = 1;
    end
    si = remove_single(si);  % Just remove all the isolated single '0' pixels
    SIS{ii} = si;
    area_sum = area_sum + sum(sum(si));
end

clear x y si tc theta ii jj
% return

%% add the clusters onto the centers
% -------------------------------------------------------------------------
% Same order:
% cl:       list of the centers' coordinates
% nd:       list of nearest distance
% Same order:
% areas:    list of all areas
% SIS:      list of all realizations of single clusters
%
% From now on, the order of cl and SIS should be coupled
% -------------------------------------------------------------------------
% Generate a temp image with buff region
temp_img = zeros(L+L/ps,L+L/ps);

% Generate an initial image
disp('Now generating an initial image...')
for ii = 1:1:N
    x = cl(ii,1);
    y = cl(ii,2);
    temp_img( x:x+L/ps-1 , y:y+L/ps-1 ) = temp_img( x:x+L/ps-1 , y:y+L/ps-1 ) + SIS{ii};
end

% %% Adjust the image using SA to eliminate overlap
% disp('Now eliminating the overlap between clusters ...')
% T = 100;
% Pr = 0.5;
% % cVF = sum(sum(temp_img1))/L/L;  % cVF: current volumn fraction
% cOL = length( find( temp_img > 1 ) );  % cOL: currenle overlap
% 
% cOL_old = cOL;
% % blk_edg_num_old = blk_edg_num;
% while T>0
%     
%     disp(T)
%     for ii = 1:1:N/2
%         
%         temp_img_old = temp_img;
%         
%         % find a cluster having overlap
%         no_list = randperm(N);
%         for jj = 1:1:N
%             pick_no_1 = no_list(ii);
%             x1 = cl(pick_no_1,1);
%             y1 = cl(pick_no_1,2);
%             if ~isempty( find( temp_img( x:x+L/ps-1 , y:y+L/ps-1 ) > 1, 1 ) ) % have overlap 
%                 break;
%             end
%         end
%         
%         pick_no_2 = ceil( rand()*N );
%         x2 = cl(pick_no_2,1);
%         y2 = cl(pick_no_2,2);
%         % Switch the location of the two clusters
%         temp_img( x1:x1+L/ps-1 , y1:y1+L/ps-1 ) = temp_img( x1:x1+L/ps-1 , y1:y1+L/ps-1 ) - SIS{pick_no_1};
%         temp_img( x2:x2+L/ps-1 , y2:y2+L/ps-1 ) = temp_img( x2:x2+L/ps-1 , y2:y2+L/ps-1 ) - SIS{pick_no_2};
%         temp_img( x1:x1+L/ps-1 , y1:y1+L/ps-1 ) = temp_img( x1:x1+L/ps-1 , y1:y1+L/ps-1 ) + SIS{pick_no_2};
%         temp_img( x2:x2+L/ps-1 , y2:y2+L/ps-1 ) = temp_img( x2:x2+L/ps-1 , y2:y2+L/ps-1 ) + SIS{pick_no_1};
%         
%         cOL = length( find( temp_img > 1 ) );  % Min overlap by Max number of '1' pixels
% 
%         if cOL < cOL_old % if better (less overlap, more "1" pixels), accept
%             
%             cl(pick_no_1,1) = x2;
%             cl(pick_no_1,2) = y2;
%             cl(pick_no_2,1) = x1;
%             cl(pick_no_2,2) = y1;
%             
% %             tswitch = SIS{pick_no_1};
% %             SIS{pick_no_1} = SIS{pick_no_2};
% %             SIS{pick_no_2} = tswitch;
%             
%             cOL_old = cOL;
%             
% %             tswitch = img_para( pick_no_1 , 3:5 );
% %             img_para( pick_no_1 , 3:5 ) = img_para( pick_no_2 , 3:5 );
% %             img_para( pick_no_2 , 3:5 ) = tswitch;
%             
% %             tswitch = areas(pick_no_1);
% %             areas(pick_no_1) = areas(pick_no_2);
% %             areas(pick_no_2) = tswitch;
%             
% %             tswitch = comps(pick_no_1);
% %             comps(pick_no_1) = comps(pick_no_2);
% %             comps(pick_no_2) = tswitch;
%             
%         else  % if worse, reject
%                 temp_img = temp_img_old;
%         end
%         
%     end
%     
%     T = T - 1;
%     Pr = Pr - (200-T)*Pr/200;
%     
%     VF_old = sum(sum( ceil(temp_img_old/(N+1)) ))/L/L;
%     VF_new = sum(sum( ceil(temp_img/(N+1)) ))/L/L;
%     disp(VF_new);
%     if abs(VF_old - VF_new) < 0.001
%         break;
%     end
%      
% end
% 
% img_para(:,1:2) = cl;
% clear tswitch

% image is the pixelated reconstruction (without adjust VF)
image = temp_img(L/ps/2+1:L+L/ps/2 , L/ps/2+1:L+L/ps/2);
image = ceil( image/(N+1) );

%% VF adjustment: radius perturbation for smooth-profile reconstruction
% Adjust the VF of smooth-profile reconstruction
vf_recon = sum( image(:) )/L^2;

img_para( :, 3:4 ) = img_para( :, 3:4 ) / sqrt( vf_recon / VF );

% if abs( vf_recon - VF ) > 0.001
% 
%     if vf_recon < VF
% 
%         area_add = L^2 * VF - sum( image(:) );
%         area_mean = area_add/N;
%         for ii = 1:1:N
%             areas(ii) = area(ii) + area_mean;
%             A = areas(ii);
%             elong = comps(ii);
%             a = sqrt( A/pi/elong );
%             b = a*elong;
%             
%             img_para(ii,3) = a;
%             img_para(ii,4) = b;
%             
%         end
%         
%     else  % vf_recon > VF
% 
%         area_substract = sum( image(:) ) - L^2 * VF;
%         area_mean = area_substract/N;
%         
%         for ii = 1:1:N
%             
%             if areas(ii) > area_mean
%                 areas(ii) = area(ii) + area_mean;
%                 A = areas(ii);
%                 elong = comps(ii);
%                 a = sqrt( A/pi/elong );
%                 b = a*elong;
% 
%                 img_para(ii,3) = a;
%                 img_para(ii,4) = b;       
%             end
%             
%         end
%         
%         
%     end
% 
% end

%% VF adjustment: pixel adding/substracting for pixelated reconstruction
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Use following part for adjust VF of the pixelated reconstruction
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %% Final tune: add/substract pixels to cluster edge to satisfy the VF
% AN = VF*L^2 - sum(sum(image));  % AN: add number
% if AN > 0
%     C = find_edge_0(image);
%     no_list = randperm(AN);
%     for ii = 1:1:min(AN,length(C))
%         no = no_list(ii);
%         x = C(no,1);
%         y = C(no,2);
%         image(x,y) = 1;
%     end
% else
%     AN = abs(AN);
%     C = find_edge_1(image);
%     no_list = randperm(AN);
%     for ii = 1:1:min(AN,length(C))
%         no = no_list(ii);
%         x = C(no,1);
%         y = C(no,2);
%         image(x,y) = 0;
%     end
% end
% -------------------------------------------------------------------------

