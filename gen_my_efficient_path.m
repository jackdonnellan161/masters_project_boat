%% John Donnellan
function [waypoints,cost] = gen_my_efficient_path(field,boat)
        dx = diff(field.xs(1:2));
        dy = diff(field.ys(1:2));
        E3 = E3_costs(field.U,field.V,dx,dy);
        xy = combvec(field.xs,field.ys)';
        s = [find(field.xs == boat.x0(1)),find(field.ys == boat.x0(2))]; 
        e = [find(field.xs == boat.xd(1)),find(field.ys == boat.xd(2))];
        startpoint = sub2ind([length(field.xs) length(field.ys)],s(1),s(2));
        endpoint = sub2ind([length(field.xs) length(field.ys)],e(1),e(2));
        [cost,path] = dijkstra(xy,E3,startpoint,endpoint);
        waypoints = xy(path,:);
end