function [Isec_out,Isec_valid] = Subfun_GetLineIntersection(ca, polygon, p0, p1, i, q0, q1, j,rn)
%% calculate the intersection of line segments (p0,p1) and (q0,q1) and judge whether the intersection is valid to form a splitting line with p1
    % ca: ca=2 represents the extended line of (p0,p1) intersects with a side (vertices not included) of the polygon; ca=2 represents the extended line of (p0,p1) intersects with a vertex of the polygon
    % polygon: the vertices of the the polygon
    % p0: the point before the concave point
    % p1: concave point
    % i: indicate the location of the concave point
    % q0: the start point of a intersecting line segment
    % q1: the end point of a intersecting line segment
    % j: indicate the location of the intersecting line segment
    % rn: the number of decimal places to be retained
    % Isec_out: the coordinate of the intersection
    % Isec_valid: judge whether the split line is in the polygon (1, yes; 0, no)
%%
    n = size(polygon,1);
    % caculate the coordinate of the intersection
    if ca == 1
        if (p0(1) == p1(1)) 
            k2 = (q1(2) - q0(2)) / (q1(1) - q0(1));
            b2 = q0(2) - k2 * q0(1);
            Isec = [p0(1), k2 * p0(1) + b2]; 
        elseif(q0(1) == q1(1))
            k1 = (p1(2) - p0(2)) / (p1(1) - p0(1));
            b1 = p0(2) - k1 * p0(1);
            Isec = [q0(1), k1 * q0(1) + b1];
        else 
            k1 = (p1(2) - p0(2)) / (p1(1) - p0(1));
            b1 = p0(2) - k1 * p0(1);
            k2 = (q1(2) - q0(2)) / (q1(1) - q0(1));
            b2 = q0(2) - k2 * q0(1);
            x = (b2 - b1) / (k1 - k2);
            y = k1 * x + b1;
            Isec = [x, y];
        end

        % find a approximate coordinate with "rn" decimal places and make the polyline (p0,p1,Isec) non-concave
        find_Isec = 0;
        vec1 = p1-p0;
        Isec_temp = [floor(Isec*10^rn);round(Isec*10^rn);ceil(Isec*10^rn)]/10^rn;
        for loop1 = 1:3
            for loop2 = 1:3
                vec2 = [Isec_temp(loop1,1),Isec_temp(loop2,2)]-p1;
                if vec1(1)*vec2(2) - vec1(2)*vec2(1) >= 0
                    find_Isec = 1;
                    break;
                end
            end
            if (find_Isec ==1)
                break;
            end
        end
        Isec_out = [Isec_temp(loop1,1),Isec_temp(loop2,2)];

        %judge whether the split line in the polygon
        Isec_valid = 1;
        for loop2 = 0:n-1
            if (loop2==i || loop2==i+1 || loop2==j)
                continue;
            else
                lineseg_A  = [polygon(loop2+1, 1) polygon(loop2+1, 2);polygon(mod(loop2+1,n)+1, 1) polygon(mod(loop2+1,n)+1, 2)];
                lineseg_B  = [p1(1) p1(2);Isec_out(1) Isec_out(2)]; 
                [xi,~] =polyxpoly(lineseg_A(:,1),lineseg_A(:,2),lineseg_B(:,1),lineseg_B(:,2));
                if ~isempty(xi)
                    Isec_valid = 0;
                    Isec_out = [];
                    break;
                end
            end
        end
    elseif ca==2
        Isec_out = q0;

        %judge whether the split line in the polygon
        Isec_valid = 1;
        for loop2 = 0:n-1
            if (loop2==i || loop2==i+1 || loop2==j-1 || loop2==j)
                continue;
            else
                lineseg_A  = [polygon(loop2+1, 1) polygon(loop2+1, 2);polygon(mod(loop2+1,n)+1, 1) polygon(mod(loop2+1,n)+1, 2)]; 
                lineseg_B  = [p1(1) p1(2);Isec_out(1) Isec_out(2)]; 
                [xi,~] =polyxpoly(lineseg_A(:,1),lineseg_A(:,2),lineseg_B(:,1),lineseg_B(:,2));
                if ~isempty(xi)
                    Isec_valid = 0;
                    Isec_out = [];
                    break;
                end
            end
        end
    end


end
