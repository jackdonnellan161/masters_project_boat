%% John Donnellan
%% Run Monte-Carlo Runs with random start and end positions in a gyre
close all; clear all; clc

rng(100)
%% Init Field
field.xs = [0:0.01:2]; %m
field.ys = [0:0.01:1]; %m
field.rho_water = 1000;
field.A_field = 0.01;
field.a = 0;
field.b = 1;

%Sizes are correct, think row/column vs x/y
[field.U,field.V,loadables.U_field,loadables.V_field] = velocity_field_init('gyre', field);

%% Init Boat
boat = boat_init('asym_point_particle',field.rho_water);

%% Init Control System
[controls, actuators, nav] = control_system_init('PI', 'direct_thrust', 'nav_passthrough');

%% Init Sim
%Sim Parameters
simu.dt = 0.005;
simu.xdot_d = 0.10;
simu.allow_disable = 0;
simu.thrust_disable_thresh = 0; 
simu.rot_disable_thresh = 0.0;
simu.cte_disable_thresh = 100;

%Constants
constants.rtod = 180/pi;
constants.dtor = pi/180;


for ii = 1:100
    clc
    % Get Start and End Points
    [start_points, end_points] = monte_carlo_get_points(1,field)
    
    %Boat IC Definition
    % boat.x0 = waypoints_load(1,:); %m
    boat.x0 = start_points(1,:); %m, overwritten if 'user_defined' waypoints
    boat.xd = end_points(1,:); %m, overwritten if 'user_defined' waypoints
    boat.v0 = [0 0]; %m/s
    boat.v0_mag = norm(boat.v0);
    boat.a0 = [0 0]; %m/s^2
    boat.alpha0 = 0;
    boat.theta0 = -90.*pi/180; %rad, init heading
    boat.theta_dot0 = 0; %rad/s
    boat.theta_ddot0 = 0; %rad/s^2

    % Waypoints
    %     waypoints_str = waypoint_str;
    %     switch waypoints_str
    %         case 'user_defined'
    %             % Put initial waypoint close enough that boat won't try to go
    %             % to it, but not so close that we get a divide by zero error in
    %             % the cte calc (fix this and remove when possible)
    %             wp_offset = [diff(field.xs(1:2)) diff(field.ys(1:2))];
    %             if (boat.x0(1) == max(field.xs))
    %                 wp_offset(1) = -diff(field.xs(1:2));
    %             end
    %             if (boat.x0(2) == max(field.ys))
    %                 wp_offset(2) = -diff(field.ys(1:2));
    %             end
    %             simu.waypoints = [boat.x0 + wp_offset; boat.xd];
    %         case 'my_efficient_path'
    %             [simu.waypoints,simu.predicted_cost] = gen_my_efficient_path(field,boat);
    %         case 'ref_efficient_path'
    % 
    %     end

    % Run straight line paths
    wp_offset = [diff(field.xs(1:2)) diff(field.ys(1:2))];
    if (boat.x0(1) == max(field.xs))
        wp_offset(1) = -diff(field.xs(1:2));
    end
    if (boat.x0(2) == max(field.ys))
        wp_offset(2) = -diff(field.ys(1:2));
    end
    simu.waypoints = [boat.x0 + wp_offset; boat.xd];
    sim('boat_sim_git.slx')
    simout_straight(ii) = simout;
    T_tot_straight(ii) = simout_straight(ii).max();
    
    %Run calculated efficient paths
    [simu.waypoints,simu.predicted_cost] = gen_my_efficient_path(field,boat);
    sim('boat_sim_git.slx')
    simout_efficient(ii) = simout;
    T_tot_efficient(ii) = simout_efficient(ii).max();
end

T_tots = [T_tot_straight', T_tot_efficient'];