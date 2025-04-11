function [splitted_poly,splitted_poledges,issplitted_out] = Subfun_SplitConcavePolygon(polygon,poledges,rn)
%% Split concave polygons by using vector method
    % polygon: the vertices of a sub-polygon, which is a part of the original entire polygon
    % poledges: the inner/outer edge labels of the sub-polygon
    % rn: the number of decimal places to be retained
    % splitted_poly: the splitted polygons of the sub-polygon
    % splitted_poledges: the inner/outer edge labels of the the sub-polygon in regard to the original entire polygon
    % issplitted: issplitted = 1 represents that the sub-polygon has been
    % splitted to 2 new sub-polygons; issplitted = 0 represents that the sub-polygon is already convex
%%
    n = size(polygon, 1);
    issplitted = 0;
    subsplitted_flag = 0;
    for i=0:n-1
        p0 = polygon(i+1, :); %the point before the checking point
        p1 = polygon(mod(i+1,n)+1, :); % the checking point
        p2 = polygon(mod(i+2,n)+1, :); % the point after the checking point
        vec1 = p1 - p0;
        vec2 = p2 - p1;
        %judge the checking point whether it is a concave point
        if (vec1(1)*vec2(2) - vec1(2)*vec2(1) < 0)%it is a concave point
            for j=0:n-1
                if (j == i-1 || j == i || j == i + 1) %skip segments (q0,q1) that cannot intersect with the extended line of segment (p0,p1)
                    continue;
                end
                q0 = polygon(j+1, :);%the start point of a line segment under search
                q1 = polygon(mod(j+1,n)+1, :); %the end point of a line segment under search
                
                vec3 = q0 - p0;
                vec4 = q1 - p0;
                % judge whether the extended line of segment (p0,p1) intersects with segments (q0,q1)
                if ((vec1(1)*vec3(2) - vec1(2)*vec3(1)) * (vec1(1)*vec4(2) - vec1(2)*vec4(1)) < 0)%the extended line of segment (p0,p1) intersects with segments (q0,q1)
                    [intersection,Isec_valid] = Subfun_GetLineIntersection(1,polygon, p0, p1,i, q0, q1, j,rn);%calculate the intersection and judge whether the intersection is valid to form a splitting line with p1
                    if Isec_valid==1
                        if (i < j)
                            %split to 2 polygons
                            polygon1 = [polygon(i+2:j+1,:);intersection];
                            polygon1_poledges = [poledges(i+2:j+1),1]; % mark the inner edge
                            polygon2 = [polygon(1:i+2,:);intersection;polygon(j+2:n,:)];
                            polygon2_poledges = [poledges(1:i+1),1,poledges(j+1:n)];
                            
%                             polygon1 = [polygon(i+2:j+1,:);intersection];
%                             polygon1_poledges = [poledges(i+2:j+1),1]; % mark the inner edge
%                             polygon2 = [polygon(i+2,:);intersection;polygon(j+2:n,:);polygon(1:i+1,:)];
%                             polygon2_poledges = [1,poledges(j+1:n),poledges(1:i+1)];
                        else
                            polygon1 = [intersection;polygon(j+2:i+2,:)];
                            polygon1_poledges = [poledges(j+1:i+1),1];
                            polygon2 = [polygon(1:j+1,:);intersection;polygon(i+2:n,:)];
                            polygon2_poledges = [poledges(1:j+1),1,poledges(i+2:n)];
                            
%                             polygon1 = [polygon(i+2:n,:);polygon(1:j+1,:);intersection];
%                             polygon1_poledges = [poledges(i+2:n),poledges(1:j+1),1];
%                             polygon2 = [polygon(i+2,:);intersection;polygon(j+2:i+1,:)];
%                             polygon2_poledges = [1,poledges(j+1:i+1)];
                        end
                    splitted_poly = {polygon1,polygon2};
                    splitted_poledges = {polygon1_poledges,polygon2_poledges};
                    issplitted = 1;
                    break;
                    else
                        continue;
                    end

                %special case£¬the extended line intersects with some certain vertex q0 of the polygon
                elseif((vec1(1)*vec3(2) - vec1(2)*vec3(1)) == 0)
                    [~,Isec_valid] = Subfun_GetLineIntersection(2,polygon, p0, p1, i, q0, q1, j,rn);%judge the effectiveness of the intersection
                    if Isec_valid==1
                        if (i < j)
                            polygon1 = polygon(i+2:j+1,:);
                            polygon1_poledges = [poledges(i+2:j),1]; % mark the inner edge
                            polygon2 = [polygon(1:i+2,:);polygon(j+1:n,:)];
                            polygon2_poledges = [poledges(1:i+1),1,poledges(j+1:n)];

%                             polygon1 = polygon(i+2:j+1,:);
%                             polygon1_poledges = [poledges(i+2:j),1]; % mark the inner edge
%                             polygon2 = [polygon(i+2,:);polygon(j+1:n,:);polygon(1:i+1,:)];
%                             polygon2_poledges = [1,poledges(j+1:n),poledges(1:i+1)];
                        else
                            polygon1 = polygon(j+1:i+2,:);
                            polygon1_poledges = [poledges(j+1:i+1),1];
                            polygon2 = [polygon(1:j+1,:);polygon(i+2:n,:)];
                            polygon2_poledges = [poledges(1:j),1,poledges(i+2:n)];

%                             polygon1 = [polygon(i+2:n,:);polygon(1:j+1,:)];
%                             polygon1_poledges = [poledges(i+2:n),poledges(1:j),1];
%                             polygon2 = [polygon(i+2,:);polygon(j+1:i+1,:)];
%                             polygon2_poledges = [1,poledges(j+1:i+1)];
                        end
                        splitted_poly = {polygon1,polygon2};
                        splitted_poledges = {polygon1_poledges,polygon2_poledges};
                        issplitted = 1;
                        break
                    else
                        continue;
                    end
                end      
            end
        end
        if issplitted ==1
            break;%found a effective splitting point
        end
    end
    
    %secondly split£¬to meet Geometry.f90 polygon_contains_point_2d_convex ( n, v, p,inside )
    if issplitted ==0%if there's no effective concave splitting point found
        [subsplitted_poly,subsplitted_edges,subsplitted_flag] = Subfun_CheckRequirementofGeometryf90(polygon,poledges);
        if (subsplitted_flag==1)
            splitted_poly = subsplitted_poly;
            splitted_poledges = subsplitted_edges;
        else
            splitted_poly = {polygon};
            splitted_poledges = {poledges};
        end
    end
    
    issplitted_out = (issplitted | subsplitted_flag);
end

