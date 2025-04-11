function [subsplitted_poly,subsplitted_edges,subsplitted_flag] = Subfun_CheckRequirementofGeometryf90(splitted_poly,poledges)
%% judge if the splitted polygons satisfy the condition in Geometry.f90 polygon_contains_point_2d_convex ( n, v, p,inside )
% further split the splitted polygon if there are more than three continuous vertices on the same line
% splitted_poly: the vertices data of the splitted polygon
% poledges: the inner/outer edge labels of the splitted polygon
% subsplitted_poly: the output vertices data after a sub-split of the input splitted polygon
% subsplitted_edges: the inner/outer edge labels of the sub-split polygons
% subsplitted_flag: 0 (not need to split); 1 (sub-splitted)
    
    singular_flag = 0;
    polycheck = splitted_poly;
    n = length(polycheck);
    subsplitted_poly = splitted_poly; % default
    subsplitted_edges = poledges; % default
    subsplitted_flag = 0;
    for loop = 3:n+2
        p0 = polycheck(loop-2,:);
        p1 = polycheck(mod(loop-2,n)+1,:);
        p2 = polycheck(mod(loop-1,n)+1,:); %这里进行切割
        p3 = polycheck(mod(loop,n)+1,:);
        k1 = (p1(2)-p0(2))/(p1(1)-p0(1));
        k2 = (p2(2)-p1(2))/(p2(1)-p1(1));
        k3 = (p3(2)-p2(2))/(p3(1)-p2(1));
        if (k2 == k3 && k1 ~= k2)%存在连续三个点共线，可能过不了Geometry.f90 polygon_contains_point_2d_convex的要求
            splitted_poly1 = [p0;p1;p2];
            splitted_poledges1 = [poledges(loop-2),poledges(mod(loop-2,n)+1),1];
            if (mod(loop-1,n)+1 > loop-2)
                splitted_poly2 = [polycheck(mod(loop-1,n)+1:n,:);polycheck(1:loop-2,:)];
                splitted_poledges2 =[poledges(mod(loop-1,n)+1:n),poledges(1:loop-3),1];
            else
                splitted_poly2 = [polycheck(mod(loop-1,n)+1:loop-2,:)];
                splitted_poledges2 =[poledges(mod(loop-1,n)+1:loop-3),1];
            end
            subsplitted_poly = {splitted_poly1,splitted_poly2};
            subsplitted_edges = {splitted_poledges1,splitted_poledges2};
            subsplitted_flag = 1;
            break;
        end
    end
end