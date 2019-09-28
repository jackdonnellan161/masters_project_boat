%Jack Donnellan
%Sim Init v1.0; 1/21/19
%Run this script to initialize simulink params
% close all; clear all; clc
xs = [0:.01:2]; %m
ys = [0:.01:1]; %m

%% Velocity Field Definition
rho_water = 1000;
A_field = 0.01;
a = 0;
b = 1;
f = a.*xs.^2 + b.*xs;
dfdx = (2.*a.*xs + b)';
U_field = (-pi.*A_field.*sin(pi.*f)'*cos(pi.*ys))'; %m/s
V_field = (pi.*A_field.*dfdx.*cos(pi.*f)'*sin(pi.*ys))'; %m/s
U_field_load = flipud(U_field);
V_field_load = flipud(V_field);
% 
% V_field= zeros(size(U_field));
% U_field= -0.1*ones(size(V_field));
% U_field_load = flipud(U_field);
% V_field_load = flipud(V_field);

% Uncomment to zero out field
U_field= zeros(size(U_field));
V_field= zeros(size(V_field));
U_field_load= zeros(size(U_field));
V_field_load= zeros(size(V_field));
quiver(xs,ys,U_field,V_field)

nx = 200;
ny = nx;
field_scale = 0.01;
b_field = 6249*field_scale;
lambda_field = 10000*field_scale;
D = .200;
F = 100;
R = 0.02;
f = ys.*10^-13;
alpha = D/R*10^-13;
gamma_field = F*pi/R/b_field;
n = 0.01;


xs = linspace(0,lambda_field,nx);
ys = linspace(0,b_field,ny);
% Psi_in = gamma_field.*(b_field/pi)^2.*sin(pi.*ys'./b_field)*...
%     [exp((xs-lambda_field).*pi/b_field) + exp(-xs.*pi./b_field)- 1];
% A = -alpha/2 + sqrt(alpha^2/4 + n^2);
% B = -alpha/2 - sqrt(alpha^2/4 + n^2);

%Field Scale = 0.01
A = 0.0001;
B = -0.2;

%Field Scale = 0.1
% A = 0.001;
% B = -0.01;

p = (1 - exp(B.*lambda_field))./(exp(A.*lambda_field)-exp(B.*lambda_field));
q = 1-p;
Psi_in = gamma_field.*(b_field/pi)^2.*sin(pi.*ys'./b_field)*...
    (p.*exp(A.*xs) + q.*exp(B.*xs)- 1);

figure
contour(xs,ys,Psi_in,10)

U_field = zeros(size(Psi_in));
V_field = zeros(size(Psi_in));


for ii = 2:length(xs)-1
    for jj = 2:length(ys)-1
        V_field(ii ,jj) = -(Psi_in (ii,jj+1) - Psi_in (ii,jj-1)) ./ (2*b_field/ny);
        U_field(ii, jj) = (Psi_in (ii+1,jj) - Psi_in (ii-1,jj)) ./ (2*lambda_field/nx);        
    end
end

%For Path Planning Script
uv = zeros(length(U_field(:,1))*length(U_field(1,:)),2);
qq = 1;
for ii = 1:length(U_field(:,1))
    for kk = 1:length(U_field(1,:))
        uv(qq,:) = [U_field(ii,kk) V_field(ii,kk)];
        qq = qq+1;
    end
end
v_max = 0.1;
v_scale_flow = v_max/max(max(max(U_field)),max(max(V_field)));
% U_field_load = v_scale_flow.*flipud(U_field);
% V_field_load = v_scale_flow.*flipud(V_field);
U_field_load = v_scale_flow.*U_field;
V_field_load = v_scale_flow.*V_field;

%% MC Setup
E3 = E3_test(U_field,V_field);
monte_carlo_get_points()
xy = combvec([1:1:length(xs)],[1:1:length(ys)])';

% figure
% contour(xs,ys,U_field)
% figure
% contour(xs,ys,V_field)
% figure
% quiver(xs,ys,U_field_load,V_field_load)

% [X,Y] = meshgrid(xs,ys);
% STARTXY = combvec(xs,ys);
% figure
% streamline(X,Y,U_field_load,V_field_load,STARTXY(1,1:2:end),STARTXY(1,1:2:end))

% quiver(xs,ys,U_field,V_field);
%% Boat Params
m_boat = 0.120; %kg
l_boat = 0.16; %m
w_boat = 0.08; %m
depth_boat = 0.06; %m
ma_x = pi*rho_water*(w_boat/2)^2*(depth_boat); %kg
ma_y = pi*rho_water*(l_boat/2)^2*(depth_boat); %kg
Izz_a = rho_water*((l_boat/2)^2-(w_boat/2)^2)^2*depth_boat; %kg * m^2
% Izz_boat = 3e-2; %kg*m
Izz_boat = 1; %kg*m
d_motors = w_boat-0.02; %distance between the centers of motor shafts, m
d_cg = [0 0 0]; %distance from IMU to CG
ax = w_boat*depth_boat; %m^2
ay = depth_boat*l_boat; %m^2
cdx = 0.1;
cdy = 0.5;
% k_drag_para = 0.3; 
% k_drag_perp = 0.7; 
k_drag_para = 1; 
k_drag_perp = 1; 
k_drag_rot = 0.4;
d_com_aero = 0.1;

%% Actuator Params (Not currently used)
L_motor = 0.1;
kphi_motor = 0.3;
J_motor = 0.1;
b_motor = 0.01;
R_motor = 2.0;

%% Control Params
kpx = 1;
kpy = 1;
kix = 0.1;
kiy = 0.1;

%% Sim Params
dt = 0.01;
use_imu_model = 0;
xd = 1.5;
yd = 0.5;

xdot_d = 0.10;
% waypoints = [1.5 0.5; 1.5 0.5];
%waypoints = [x0_boat; 1.5 0.5;1.0 0.75; 0.5 0.5;1.0 0.25];
% waypoints = [0.5 0.5; 0.7000 0.6400; 1.0000 0.3000; 1.1200 0; 1.4600 0.1200; 1.7000 0.4600;1.5 0.5];
% waypoints = waypoints./100;
% waypoints = [x0_boat;1.5 0.5];
close all
waypoints_load = [key_points(1,:); key_points(2:end,:)];
% waypoints_load = [key_points(1,:) + 0.01;key_points(end,:)];

allow_disable = 1;
% thrust_disable_thresh_upper = 0.05;
% thrust_disable_thresh_lower = 0.02;
thrust_disable_thresh = 1; 
rot_disable_thresh = 0.0;
cte_disable_thresh = 100;

%% Boat IC Definition
% x0_boat = [1 0.5]; %m
x0_boat = waypoints_load(1,:); %m
v0_boat = [0 0]; %m/s
v0_boat_mag = norm(v0_boat);
alpha0_boat = 0;
a0_boat = [0 0]; %m/s^2
theta0_boat = 0.*pi/180; %rad
theta_dot0_boat = 0; %rad/s
theta_ddot0_boat = 0; %rad/s^2

%% Airfoil Params
%alpha    CL     CD
%NACA Airfoil, reference "Aerodynamic characteristics of seven symmetrical
%airfoil sections..." by Sheldahl et al.

airfoil_lookup = ...
[0.0000 0.0000 0.0091;
1.0000 0.1100 0.0092;
2.0000 0.2200 0.0094;
3.0000 0.3300 0.0098;
4.0000 0.4400 0.0105;
5.0000 0.5500 0.0114;
6.0000 0.6600 0.0126;
7.0000 0.7390 0.0143;
8.0000 0.8240 0.0157;
9.0000 0.8946 0.0173;
10.0000 0.9440 0.0191;
11.0000 0.9572 0.0211;
12.0000 0.9285 0.0233;
13.0000 0.8562 0.0257;
14.0000 0.7483 0.0283;
15.0000 0.6350 0.0312;
16.0000 0.5384 0.1240;
17.0000 0.4851 0.2170;
18.0000 0.4782 0.2380;
19.0000 0.4908 0.2600;
20.0000 0.5247 0.2820;
21.0000 .5616 0.3050;
22.0000 0.6045 0.3290;
23.0000 0.6528 0.3540;
24.0000 0.7015 0.3790;
25.0000 0.7511 0.4050;
26.0000 0.8055 0.4320;
27.0000 0.8788 0.4600;
30.0000 0.8550 0.5700;
35.0000 0.9800 0.7450;
40.0000 1.0350 0.9200;
45.0000 1.0500 1.0750;
50.0000 01.0200 1.2150;
55.0000 0.9550 1.3450;
60.0000 0.8750 1.4700;
65.0000 0.7600 1.5750;
70.0000 0.6300 1.6650;
75.0000 0.5000 1.7350;
80.0000 0.3650 1.7800;
85.0000 0.2300 1.8000;
90.0000 0.0900 1.8000;
95.0000 -0.0500 1.7800;
100.0000 -.1850 1.7500;
105.0000 -.3200 1.7000;
110.0000 -.4500 1.6350;
115.0000 -.5750 1.5550;
120.0000 -0.6700 1.4650;
125.0000 -0.7600 1.3500;
130.0000 -0.8500 1.2250;
135.0000 -0.9300 1.0850;
140.0000 -0.9800 0.9250;
145.0000 -0.9000 0.7550;
150.0000 -0.7700 0.5750;
155.0000 -0.6700 0.4200;
160.0000 -0.6350 0.3200;
165.0000 -0.6800 0.2300;
170.0000 -0.8500 0.1400;
175.0000 -0.6600 0.0550;
180.0000 0.0000 0.0250];

aifoil_negative_side = [-airfoil_lookup(2:end,1),...
    -airfoil_lookup(2:end,2),airfoil_lookup(2:end,3)];

airfoil_lookup = [flip(aifoil_negative_side,1);airfoil_lookup];
alpha_lookup = airfoil_lookup(:,1);
cl_lookup = airfoil_lookup(:,2);
cd_lookup = airfoil_lookup(:,3);

% figure
% scatter(alpha_lookup,[cl_lookup])
% 
% figure
% scatter(alpha_lookup,[cd_lookup])

%% Constants
rtod = 180/pi;
dtor = pi/180;

close all