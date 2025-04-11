function k = Fun_Curvature(x, y)
%% Used for calculating the curvature of the input curve
    % (x,y): the input point coordinates of the curve with 2 columns.
    % k: the calcultated curvature
%%
    % Expand the coordinates to include the start and end loops.
    x_ext = [x(end); x; x(1)];
    y_ext = [y(end); y; y(1)];
    
    % First-order derivative (using central difference)
    dx = 0.5 * (x_ext(3:end) - x_ext(1:end-2));
    dy = 0.5 * (y_ext(3:end) - y_ext(1:end-2));
    
    % Second-order derivative (using central difference)
    ddx = 0.5 * (x_ext(3:end) - 2*x_ext(2:end-1) + x_ext(1:end-2));
    ddy = 0.5 * (y_ext(3:end) - 2*y_ext(2:end-1) + y_ext(1:end-2));
    
    % Calculating the curvature
    numerator = abs(dx .* ddy - dy .* ddx);
    denominator = (dx.^2 + dy.^2).^(3/2);
    k = numerator ./ denominator;
    
    k(isnan(k)) = 0;
end