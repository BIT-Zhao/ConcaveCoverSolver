function [splitted_poly_save,polyedges_save] = Fun_SplitConvaePolygon(poly,rn)
%% Convex split of concave polygon
    % poly: the vertex coordinates of the polygon
    % rn: The number of decimal places to be retained, to avoid rounding errors in the final.f90 file that may cause concave polygons.
    % splitted_poly_save: the splitted polygons of the input "poly"
    % polyedges_save: the inner/outer edge labels. 0 (the outter edge); 1 (the inner edge)
%%
    poly = round(poly,rn);
    splitted_poly_save = {}; % the saved splitted polygons
    polyedges_save = {}; % the saved edge types of the splitted polygons
    polygon = {poly};
    polyedges = {zeros(1,length(poly))};
    while ~isempty(polygon) 
        [splitted_poly,splitted_poledges,issplitted,] = Subfun_SplitConcavePolygon(polygon{1},polyedges{1},rn);
        
%         figure; hold on;
%         colors = hsv(length(splitted_poly));
%         for i = 1:length(splitted_poly)
%             part = splitted_poly{i};
%             fill(part(:,1), part(:,2), colors(i,:));
%         end
%         title('Split result');

        if issplitted == 0 %The current polygon is already a convex polygon.
            splitted_poly_save = [splitted_poly_save,splitted_poly];
            polyedges_save = [polyedges_save;splitted_poledges];
            polygon = polygon(2:end);
            polyedges = polyedges(2:end);
        else %Place the two polygons that have been splitted at the front of the queue to await the next step of split.
            polygon = [splitted_poly,polygon(2:end)];
            polyedges = [splitted_poledges,polyedges(2:end)];
        end
        
    end
    
end