function [max_dist, max_idx, max_dist_seg] = Fun_MaxDistance(seg_x, seg_y, approx_x, approx_y,ava_ind)
%% Calculates the error of the line segment from the point to the end point in the current ava_ind interval
    % seg_x: the x coordinates of points in the present segment
    % seg_y: the y coordinates of points in the present segment
    % approx_x: the x coordinates of the already determined vertex
    % approx_y: the y coordinates of the already determined vertex
    % ava_ind: vertices indicator which records the locations of the already determined vertices in the present segment
    % max_dist: the maximum distance from all points in [seg_x,seg_y] to line segments formed by vertices [approx_x,approx_y]
    % max_idx: the index of the maximum distance
    % max_dist_seg: the segment index identifying the maximum-distance segment
%%
    len = length(approx_x);
    for subseg = 1:len-1
        % End point of line segment  
        x1 = approx_x(subseg);
        y1 = approx_y(subseg);
        x2 = approx_x(subseg+1);
        y2 = approx_y(subseg+1);
        x = seg_x;
        y = seg_y;

        cross = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);
        d2 = (x2 - x1)^2 + (y2 - y1)^2;
        r = cross / d2;
        px = x1 + (x2 - x1) * r;
        py = y1 + (y2 - y1) * r;

        distances = sqrt((x - px).^2 + (py - y).^2);
        distances(cross <= 0) = sqrt((x(cross <= 0) - x1).^2 + (y(cross <= 0) - y1).^2);
        distances(cross >= d2) = sqrt((x(cross >= d2) - x2).^2 + (y(cross >= d2) - y2).^2);

        Dis(ava_ind(subseg):ava_ind(subseg+1),subseg) = distances(ava_ind(subseg):ava_ind(subseg+1));
    end
        
        [max_dist_dim2, max_idx_dim2] = max(Dis,[],2);
        % Maximum distance and its index
        [max_dist, max_idx] = max(max_dist_dim2);
        % Maximum-distance-matched segment
        max_dist_seg = max_idx_dim2(max_idx);
end