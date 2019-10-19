close all; clear all; clc
%% Local Planning Graph Point Initialization
xStep = 1;
yStep = xStep;
yy = 50;
xs = 0:xStep:2*yy;
ys = 0:yStep:yy;
thin_xs = 0:2*xStep:2*yy;
thin_ys = 0:2*yStep:yy; 
xy = combvec(xs,ys)';
thin_xy = combvec(thin_xs,thin_ys)';
%% Velocity Field Definition
% %Uniform Flow
% u = ones(length(xy(:,1)),1);
% v = zeros(length(xy(:,1)),1);
% xyuv = [xy u v];

% %Double Gyre
params.A = 1;
vel_points = [xy(:,1)/yy xy(:,2)/yy];
vel = stable_doublegyre(vel_points,params);
xyuv = [xy vel];

%% Local Planning Graph Point Initialization
thin_xyuv = zeros(1,4);
%Want Global Planning on Grid Half Size
for ii = 1:length(xyuv)
    if mod(xyuv(ii,1),2) == 0 && mod(xyuv(ii,2),2) == 0
        thin_xyuv = [thin_xyuv; xyuv(ii,:)];
    end
end
thin_xyuv = thin_xyuv(2:end,:);
% figure
% quiver(thin_xyuv(:,1),thin_xyuv(:,2),thin_xyuv(:,3),thin_xyuv(:,4))
%% Connections and Costs
% thin_E3 = E3_square(thin_xs,thin_ys,thin_xyuv);
tic
E3 = E3_square(xs,ys,xyuv);

%% Global Path Finding
s = [yy/2 yy/2]; 
e = [1.5*yy yy/2];
startpoint = sub2ind([length(xs) length(ys)],s(1),s(2));
endpoint = sub2ind([length(xs) length(ys)],e(1),e(2));
[cost,path] = dijkstra(xy,E3,startpoint,endpoint);
toc
% cost
% path
global_points = xy(path,:);

%% Local Path Planning between Waypoints
% num_waypoints = 1;
% waypoints = zeros(num_waypoints,2);
% kk = 1;
% space = round(length(global_points)./(num_waypoints+1));
% for ii = 1:num_waypoints
%     waypoints(ii,:) = global_points(ii*space,:);
% end
% key_points = [s;waypoints;e];
% points2 = zeros(1,2);
% for ii = 1:length(key_points)-1
%     tic
%     waypoint1 = sub2ind([length(xs) length(ys)],key_points(ii,1),key_points(ii,2));
%     waypoint2 = sub2ind([length(xs) length(ys)],key_points(ii+1,1),key_points(ii+1,2));
%     [cost,path] = dijkstra(xy,E3,waypoint1,waypoint2);
%     points2 = [points2; xy(path,:)];
%     toc
% end
% points2 = points2(2:end,:);
%% Path Plotting
figure
quiver(xyuv(:,1),xyuv(:,2),xyuv(:,3),xyuv(:,4))
hold on
scatter(global_points(:,1),global_points(:,2),'r')
% scatter(waypoints(:,1),waypoints(:,2),'b')
scatter(s(1),s(2),'g')
scatter(e(1),e(2),'g')
axis('equal')

figure
quiver(xyuv(:,1),xyuv(:,2),xyuv(:,3),xyuv(:,4))
hold on
scatter(points2(:,1),points2(:,2),'r')
axis('equal')