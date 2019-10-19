% close all; clear all; clc
%% Local Planning Graph Point Initialization
% xStep = 1;
% yStep = xStep;
% yy = 100;
% [xs,ys,xy] = grid_define(0,2*yy,0,yy,xStep);
% [thin_xs,thin_ys,thin_xy] = grid_define(0,2*yy,0,yy,2*xStep);
yy = lambda_field;
xStep = lambda_field/nx;
yStep = b_field/xStep;
% [xs,ys,xy] = grid_define(0,lambda_field,0,b_field,yStep);
xy = combvec(xs,ys)';
% [thin_xs,thin_ys,thin_xy] = grid_define(0,lambda_field,0,b_field,2*yy);
%% Velocity Field Definition
% Double Gyre
% params.A = 1;
% vel_points = [xy(:,1)/yy xy(:,2)/yy];
% vel = stable_doublegyre(vel_points,params);

vel = uv;
xyuv = [xy vel];

quiver(xyuv(:,1),xyuv(:,2),xyuv(:,3),xyuv(:,4))

% thin_vel_points = [thin_xy(:,1)/yy thin_xy(:,2)/yy];
% thin_vel = stable_doublegyre(thin_vel_points,params);
% thin_xyuv = [thin_xy thin_vel];

%% Local Planning Graph Point Initialization
% thin_xyuv = zeros(1,4);
% % Want Global Planning on Grid Half Size
% for ii = 1:length(xyuv)
%     if mod(xyuv(ii,1),2) == 0 && mod(xyuv(ii,2),2) == 0
%         thin_xyuv = [thin_xyuv; xyuv(ii,:)];
%     end
% end
% thin_xyuv = thin_xyuv(2:end,:);
% % figure
% % quiver(thin_xyuv(:,1),thin_xyuv(:,2),thin_xyuv(:,3),thin_xyuv(:,4))

%% Connections and Costs
% thin_E3 = E3_square(thin_xs,thin_ys,thin_xyuv);
tic
E3 = E3_square(xs,ys,xyuv);
toc

%% Global Path Finding
tic
%s = [round(yy/4)*2 round(yy/4)*2]; 
%e = [round(yy/4)*6 round(yy/4)*2];
s = [10 190]; 
e = [190 10];
startpoint = sub2ind([length(xs) length(ys)],s(1)/2,s(2)/2);
endpoint = sub2ind([length(xs) length(ys)],e(1)/2,e(2)/2);
[cost,path] = dijkstra(xy,E3,startpoint,endpoint);
% cost
% path
global_points = xy(path,:);
toc

% %% Local Planning Between Waypoints
num_waypoints = 10;
waypoints = zeros(num_waypoints,2);
kk = 1;
space = round(length(global_points)./(num_waypoints+1));
for ii = 1:num_waypoints
    waypoints(ii,:) = global_points(ii*space,:);
end
s = [s/200.*[lambda_field b_field]];
e = [e/200.*[lambda_field b_field]];
key_points = [s;waypoints;e];
% 
% points2 = zeros(1,2);
% for ii = 1:length(key_points)-1
%     tic
%     way_path = waypoint_planner(key_points(ii:(ii+1),:),xStep,yy,params);
%     points2 = [points2;way_path];
%     toc
% end
% points2 = points2(2:end,:);

%% Path Plotting
figure
quiver(xyuv(:,1),xyuv(:,2),xyuv(:,3),xyuv(:,4))
hold on
scatter(global_points(:,1),global_points(:,2),'r')
scatter(waypoints(:,1),waypoints(:,2),'b')
scatter(global_points(1,1),global_points(1,2),'g')
scatter(global_points(end,1),global_points(end,2),'g')
axis('equal')

% figure
% quiver(xyuv(:,1),xyuv(:,2),xyuv(:,3),xyuv(:,4))
% hold on
% scatter(points2(:,1),points2(:,2),'r')
% axis('equal')

figure
% scatter(points2(:,1),points2(:,2),'r')
hold on
scatter(global_points(:,1),global_points(:,2),'c')
% scatter(waypoints(:,1),waypoints(:,2),'b')
scatter(global_points(1,1),global_points(1,2),'g')
scatter(global_points(end,1),global_points(end,2),'g')
axis('equal')

close all