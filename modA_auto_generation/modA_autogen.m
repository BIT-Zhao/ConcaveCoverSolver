%% modA_autogen.m introduction
%Program function: 
    % Automated convex decomposition of explicitly closed polygons, generates inner/outer edge labels, and exports a toolbox-compatible 
    % Agen.f90 file to enable rapid solution of the thinnest covering problem for arbitrary polygons.
%Used solving toolbox of the thinnest covering problem£º
    % https://github.com/johngardenghi/bglcovering¡£
    % Reference: [1] E. G. Birgin, J. L. GardenghiºÍA. Laurain, ¡¶Asymptotic bounds on the optimal radius when covering a set with minimum radius identical balls¡·.

%%
clear all;
close all;
clc;

% Load polygon data of vertices_x and vertices_y. Optionally, a hand-painted figure can be imported to generate the polygon data by polygon_extract.m (its
% functions Fun_Curvature.m, Fun_MaxDistance.m, Fun_SegmentContour.m are also included). 
% !!! Vertex coordinate data are arranged in counterclockwise order for each polygon.
% load HP_polygon_parameters.mat;%to generate hand-painted Agen
load MIF_polygon_parameters.mat;%to generate approximate Minkowski_island_fractal Agen

% Regular polygon vertices generation
% nvert = 50;
% phi = 2*pi/nvert*(0:nvert).';
% phi_mid12 = (phi(1)+phi(2))/2;
% phi = phi-pi/2-phi_mid12;%Make the bottom edge parallel to the X-axis
% vertices_x = cos(phi);
% vertices_y = sin(phi);

% Coordinates processing: translate the polygon near the origin and normalization of the length along the x direction
vertices_x = vertices_x - mean(vertices_x(1:end-1),'all');
vertices_y = vertices_y - mean(vertices_y(1:end-1),'all');
scaling_ratio = (max(vertices_x,[],'all')-min(vertices_x,[],'all'));
vertices_x = vertices_x/scaling_ratio;
vertices_y = vertices_y/scaling_ratio;

% Split concave polygon
rn = 4;%The number of decimal places retained for coordinate values in the final Agen.f90 to be exported.
poly = [vertices_x(1:end-1),vertices_y(1:end-1)];
[splitted_poly_save,polyedges_save] = Fun_SplitConvaePolygon(poly,rn); % Function of split

%% Prepare splitted polygon data
modA_npols = length(splitted_poly_save);%the number of polygons

for loop = 1:modA_npols
    %the number of vertices of the polygons
    modA_nvpols(loop) = size(splitted_poly_save{loop},1) + 1;


    %put the X and Y coordinates of all the vertices into modA_putx and modA_puty, respectively
    modA_putx(loop,1:modA_nvpols(loop)) = [splitted_poly_save{loop}(:,1);splitted_poly_save{loop}(1,1)].';%X coordinates
    modA_puty(loop,1:modA_nvpols(loop)) = [splitted_poly_save{loop}(:,2);splitted_poly_save{loop}(1,2)].';%Y coordinates
    
    %the inner/outer edge mark
    modA_poledges(loop,1:modA_nvpols(loop)-1) = polyedges_save{loop};
end

% Rectangle containing A
modA_dbl = [min(modA_putx,[],'all'),min(modA_puty,[],'all')] - 0.1;
modA_dtr = [max(modA_putx,[],'all'),max(modA_puty,[],'all')] + 0.1;

% Visualize
figure; 
subplot(1,2,1); 
fill(poly(:,1), poly(:,2), 'b'); 
title('The original polygon');
subplot(1,2,2); hold on;
colors = hsv(modA_npols);
for i = 1:modA_npols
    part = splitted_poly_save{i};
    line = polyedges_save{i};
    fill(part(:,1), part(:,2), colors(i,:),'Linestyle','--');hold on;
    for j = 1:length(line)
        if line(j) == 0
            plot(modA_putx(i,j:j+1),modA_puty(i,j:j+1),'LineWidth', 2,'Color','k');hold on;
        end
    end
end
title('The splitted polygons');

figure;
colors = hsv(modA_npols);
for i = 1:modA_npols
    part = splitted_poly_save{i};
    line = polyedges_save{i};
    fill(part(:,1), part(:,2), colors(i,:),'Linestyle','--');hold on;
    for j = 1:length(line)
        if line(j) == 0
            plot(modA_putx(i,j:j+1),modA_puty(i,j:j+1),'LineWidth', 2,'Color','k');hold on;
        end
    end
end
grid on;
% title('X-axis size-normalized split polygons for ','interpreter','latex');
xlabel('X','interpreter','latex')
ylabel('Y','interpreter','latex')
xticks(-0.6:0.2:0.6);

%% generage Agen.f90
fid = fopen('Agen.f90', 'w');
fprintf(fid, 'module Agen\n  implicit none\n\n  integer :: modA_npols\n');
fprintf(fid, '  parameter ( modA_npols = %d )', modA_npols);fprintf(fid, ' !number of polygons\n\n');

fprintf(fid, '  !number of vertices (or borders) of each polygon\n');
fprintf(fid, '  integer :: modA_nvpols(modA_npols)\n');
fprintf(fid, '  data modA_nvpols(1:modA_npols) / &\n');
N_data_row = 20;
for j = 1:ceil(modA_npols/N_data_row)%Exclude the last point showing the closed polygon; fortran required 72 column characters per row, and we put 20 data per row
    if (j==ceil(modA_npols/N_data_row))            
        fprintf(fid, '       %s /\n\n',strjoin(string(modA_nvpols(1+N_data_row*(j-1):modA_npols)), ', '));
    else
        fprintf(fid, '       %s, &\n',strjoin(string(modA_nvpols(1+N_data_row*(j-1):N_data_row*j)), ', '));
    end
end
    
    
    

fprintf(fid, '  !Acquired coordinates of each boundary points\n');
fprintf(fid, '  real(kind=8) :: modA_putx(%d, %d),  modA_puty(%d, %d)\n',modA_npols,max(modA_nvpols),modA_npols,max(modA_nvpols));
N_data_row = 5;
for i = 1:modA_npols
    if (i==1)
        fprintf(fid, '  data modA_putx(%d,1:%d) / &\n',i,modA_nvpols(i));
    else
        fprintf(fid, '       modA_putx(%d,1:%d) / &\n',i,modA_nvpols(i));
    end
    for j = 1:ceil(modA_nvpols(i)/N_data_row)%fortran requires 72 column characters per row, and we put N_data_row of data per row
        if (i == modA_npols && j==ceil(modA_nvpols(i)/N_data_row))            
            elements_with_d0 = arrayfun(@(x) sprintf('%.4fd0', x), modA_putx(i, 1+N_data_row*(j-1):modA_nvpols(i)), 'UniformOutput', false);
            fprintf(fid, '       %s /\n\n',strjoin(elements_with_d0, ', '));
        elseif ( j==ceil(modA_nvpols(i)/N_data_row))            
            elements_with_d0 = arrayfun(@(x) sprintf('%.4fd0', x), modA_putx(i, 1+N_data_row*(j-1):modA_nvpols(i)), 'UniformOutput', false);
            fprintf(fid, '       %s /, &\n',strjoin(elements_with_d0, ', '));
        else
            elements_with_d0 = arrayfun(@(x) sprintf('%.4fd0', x), modA_putx(i, 1+N_data_row*(j-1):N_data_row*j), 'UniformOutput', false);
            fprintf(fid, '       %s, &\n',strjoin(elements_with_d0, ', '));
        end
    end 
end
for i = 1:modA_npols
    if (i==1)
        fprintf(fid, '  data modA_puty(%d,1:%d) / &\n',i,modA_nvpols(i));
    else
        fprintf(fid, '       modA_puty(%d,1:%d) / &\n',i,modA_nvpols(i));
    end
    for j = 1:ceil(modA_nvpols(i)/N_data_row)%fortran requires 72 column characters per row, and we put N_data_row of data per row
        if (i == modA_npols && j==ceil(modA_nvpols(i)/N_data_row))            
            elements_with_d0 = arrayfun(@(x) sprintf('%.4fd0', x), modA_puty(i, 1+N_data_row*(j-1):modA_nvpols(i)), 'UniformOutput', false);
            fprintf(fid, '       %s /\n\n',strjoin(elements_with_d0, ', '));  
        elseif ( j==ceil(modA_nvpols(i)/N_data_row))            
            elements_with_d0 = arrayfun(@(x) sprintf('%.4fd0', x), modA_puty(i, 1+N_data_row*(j-1):modA_nvpols(i)), 'UniformOutput', false);
            fprintf(fid, '       %s /, &\n',strjoin(elements_with_d0, ', '));
        else
            elements_with_d0 = arrayfun(@(x) sprintf('%.4fd0', x), modA_puty(i, 1+N_data_row*(j-1):N_data_row*j), 'UniformOutput', false);
            fprintf(fid, '       %s, &\n',strjoin(elements_with_d0, ', '));
        end
    end
end


fprintf(fid, '  !Edges ID\n');
fprintf(fid, '  integer :: modA_ploedges(%d, %d)\n',modA_npols,max(modA_nvpols)-1);
N_data_row = 20;
for i = 1:modA_npols
    if (i==1)
        fprintf(fid, '  data modA_ploedges(%d,1:%d) / &\n',i,modA_nvpols(i)-1);
    else
        fprintf(fid, '       modA_ploedges(%d,1:%d) / &\n',i,modA_nvpols(i)-1);
    end
    for j = 1:ceil((modA_nvpols(i)-1)/N_data_row)%Exclude the last point showing the closed polygon; fortran required 72 column characters per row, and we put 20 data per row
        if (i==modA_npols && j==ceil((modA_nvpols(i)-1)/N_data_row))            
            fprintf(fid, '       %s /\n\n',strjoin(string(modA_poledges(i,1+N_data_row*(j-1):modA_nvpols(i)-1)), ', '));
        elseif (j==ceil((modA_nvpols(i)-1)/N_data_row))            
            fprintf(fid, '       %s /, &\n',strjoin(string(modA_poledges(i,1+N_data_row*(j-1):modA_nvpols(i)-1)), ', '));
        else
            fprintf(fid, '       %s, &\n',strjoin(string(modA_poledges(i,1+N_data_row*(j-1):N_data_row*j)), ', '));
        end
    end
end


fprintf(fid, '  ! Rectangle containing A\n');
fprintf(fid, '  real(kind=8) :: modA_dbl(2), modA_dtr(2)\n');
fprintf(fid, '  data modA_dbl(1:2) / %.4fd0, %.4fd0 /, &\n', modA_dbl(1), modA_dbl(2));
fprintf(fid, '       modA_dtr(1:2) / %.4fd0, %.4fd0 /\n\n', modA_dtr(1), modA_dtr(2));   
        
        
fprintf(fid, 'end module Agen\n');
fclose(fid);