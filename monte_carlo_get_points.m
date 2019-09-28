num_points = 1;
for ii = 1:num_points
    x1(ii,:) = unidrnd(nx);
    y1(ii,:) = unidrnd(ny);
    
    x2(ii,:) = unidrnd(nx);
    y2(ii,:) = unidrnd(ny);
end

start_points = [x1 y1];
end_points = [x2 y2];
% 
% figure
% hold on
% xlim([0 nx])
% ylim([0 ny])
for ii = 1:num_points
   P1 = [start_points(ii,1) start_points(ii,2)];
   P2 = [end_points(ii,1) end_points(ii,2)];
   D = P2 - P1;
%    quiver( P1(1), P1(2), D(1), D(2), 0 )
end