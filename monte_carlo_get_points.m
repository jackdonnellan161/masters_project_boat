function [start_points, end_points] = monte_carlo_get_points(num_points,field)
    nx = length(field.xs);
    ny = length(field.ys);
    x1 = unidrnd(nx,num_points,1);
    y1 = unidrnd(ny,num_points,1);
    x2 = unidrnd(nx,num_points,1);
    y2 = unidrnd(ny,num_points,1);
    start_points = [field.xs(x1)' field.ys(y1)'];
    end_points = [field.xs(x2)' field.ys(y2)'];
end