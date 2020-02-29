function [start_points, end_points] = monte_carlo_get_points(num_points,field)
    
    nx = length(field.xs);
    ny = length(field.ys);
    for ii = 1:num_points
        x1(ii,:) = unidrnd(nx);
        y1(ii,:) = unidrnd(ny);
        x2(ii,:) = unidrnd(nx);
        y2(ii,:) = unidrnd(ny);
    end
    start_points = [field.xs(x1) field.ys(y1)];
    end_points = [field.xs(x2) field.ys(y2)];
end