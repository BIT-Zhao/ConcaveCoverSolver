%% polygon_extract.m
% Program function: This program is designed to extract the contour of the enclosed figure inputted and convert it into an approximate polygon.
close all;
clear all;
clc;

%% Image preprocessing and contour extraction
% Read the image and perform binarization
% img = imread('test_ori.png');
% img = imread('test_ori2.png');
% % img = imread('test_ori3.png');
% img = imread('test_ori4.png');% Read the input image.The image can be hand-painted as a simple example. Or we can directly generate in MATLAB.
% img = imread('handpainted.jpg');
img = imread('Minkowski_island_fractal.png');
grayImg = rgb2gray(img);
bwImg = imbinarize(grayImg); % Binarization

% Extract the contour points of the closed figure (in clockwise order)
boundaries = bwboundaries(bwImg);% Extract the contour of the curve
contour = boundaries{2};% Obtain the outer contour of the curve. For multiple or complex closed figures, please check the index.
x = contour(:,2);  
y = -contour(:,1); 
x = flipud(x);%in clockwise order
y = flipud(y);%in clockwise order

% Due to pixel approximation error, smoothing.
% x = smooth(x,3); % The number of smoothing points has an impact on the final polygon point count. Adjust the number of smoothing points to obtain the desired output polygon.
% y = smooth(y,3); % In our simulations, we commented it out for Minkowski_island_fractal case

%% Curvature analysis and key point detection
% Curvature calculation
k = Fun_Curvature(x(1:end-1), y(1:end-1));

[k3,ind] = maxk(k,3); % At least 3 points to generate a polygon.

% Detect high-curvature points (corner points) and low-curvature points (curve segments)
threshold = k3(3)*1.0;  % Adjusting the proportion should not exceed 1; selection of the initial point.
corner_idx = find(k >= threshold);
curve_segments = Fun_SegmentContour(corner_idx, length(x));

%% Vertices optimization
vertices_x_init = x(corner_idx);%Initializes the vertices collection
vertices_y_init = y(corner_idx);
vertices_x = [];
vertices_y = [];

% the maximum allowable error /pixels
max_error =5.0;

% Iterative optimization
for i = 1:length(curve_segments)
    seg = curve_segments{i};
    seg_x = x(seg);
    seg_y = y(seg);
    
    approx_x = [seg_x(1), seg_x(end)];% Initial approximation: Join the beginning and end of a curve segment
    approx_y = [seg_y(1), seg_y(end)];
    
    ava_ind = [1,length(seg)];%Initial vertices indicator: an interval that can be used to distribute vertices for approx_x and approx_y
    
    [current_error, max_idx, max_dist_seg] = Fun_MaxDistance(seg_x, seg_y, approx_x, approx_y,ava_ind);% Calculates the error of the line segment from the point to the end point in the current ava_ind interval

    % Increments the vertices until the error requirement is reached
    while current_error > max_error
        approx_x = [approx_x(1:max_dist_seg),seg_x(max_idx),approx_x(max_dist_seg+1:end)];% add new vertices
        approx_y = [approx_y(1:max_dist_seg),seg_y(max_idx),approx_y(max_dist_seg+1:end)];
        
        ava_ind = [ava_ind(1:max_dist_seg),max_idx,ava_ind(max_dist_seg+1:end)];% update vertices indicator, which records the locations of the already determined vertices in the present segment
        
        [current_error, max_idx, max_dist_seg] = Fun_MaxDistance(seg_x, seg_y, approx_x, approx_y,ava_ind);% update maximum distance error; newly added vertex and its loction indicator
    end
    
    % Adds the optimized vertices to the collection
    vertices_x = [vertices_x; approx_x(1:end-1).'];
    vertices_y = [vertices_y; approx_y(1:end-1).'];
end

% Explicitly closed polygon coordinates
vertices_x = [vertices_x;vertices_x(1)];
vertices_y = [vertices_y;vertices_y(1)];

%% test
% vertices_x = [vertices_x(1:7);vertices_x(7);vertices_x(1)];
% vertices_y = [vertices_y(1:7);vertices_y(end-1);vertices_y(1)];

% vertices_x = [vertices_x(1:15);vertices_x(16);vertices_x(1)];
% vertices_y = [vertices_y(1:15);vertices_y(16);vertices_y(1)];
% 
% vertices_x = [vertices_x(1:16);vertices_x(17);vertices_x(1)];
% vertices_y = [vertices_y(1:16);vertices_y(17);vertices_y(1)];

% vertices_x = [vertices_x(2:16);vertices_x(17);vertices_x(1:2)];
% vertices_y = [vertices_y(2:16);vertices_y(17);vertices_y(1:2)];



% vertices_x = [vertices_x(1:2);vertices_x(16);vertices_x(16);(vertices_x(1)+vertices_x(16))/2;(vertices_x(1)+vertices_x(16))/2;vertices_x(1)];%target
% vertices_y = [vertices_y(1:2);vertices_y(2);vertices_y(16);vertices_y(16);vertices_y(1);vertices_y(1)];

% vertices_x = [vertices_x(2);vertices_x(16);vertices_x(16);(vertices_x(1)+vertices_x(16))/2;(vertices_x(1)+vertices_x(16))/2;vertices_x(1:2)];%target
% vertices_y = [vertices_y(2);vertices_y(2);vertices_y(16);vertices_y(16);vertices_y(1);vertices_y(1:2)];

% vertices_x = [0;1;2;1;1;0;0];%target
% vertices_y = [0;0;2;2;1;1;0];

% vertices_x = [0;2;2;1;1;0;0];%target
% vertices_y = [0;0;1;1;2;2;0];

% vertices_x = [vertices_x(1:3);vertices_x(3);80;vertices_x(1)];
% vertices_y = [vertices_y(1:3);vertices_y(1);vertices_y(1);vertices_y(1)];
% 
% vertices_x = [vertices_x(1:2);vertices_x(7);vertices_x(7);vertices_x(1)];
% vertices_y = [vertices_y(1:2);vertices_y(7);vertices_y(end-1);vertices_y(1)];
%% Visualize
figure;
plot(x, y, 'b-', 'LineWidth', 1);grid on;
hold on;
plot(vertices_x, vertices_y, 'ro--', 'LineWidth', 1, 'Markersize',2.5);
leg = legend('$\mathcal{S}_{3}$', '$\tilde{\mathcal{S}}_{3}$','interpreter','latex');
leg.ItemTokenSize = [10,10];
% title('Edge of $\mathcal{S}_{3}$ and $\tilde{\mathcal{S}}_{3}$','interpreter','latex');
xlabel('X/pixel','interpreter','latex');
ylabel('Y/pixel','interpreter','latex');
%% save polygon parameters
save("MIF_polygon_parameters","vertices_x","vertices_y");
% save("HP_polygon_parameters","vertices_x","vertices_y");

