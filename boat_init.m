%% John Donnellan
%  Consolidated boat model parameters with switch for different aero models

function boat = boat_init(aero_case,rho)
%% Size/Dimensions
boat.length = 0.16; %m
boat.width = 0.08; %m
boat.depth = 0.06; %m
boat.d_motors = boat.width-0.02; %distance between the centers of motor shafts, m
boat.Ax = boat.width*boat.depth; %m^2
boat.Ay = boat.depth*boat.length; %m^2

%% Mass/Inertias
boat.m = 0.120; %kg
boat.added_mass_x = pi*rho*(boat.width/2)^2*(boat.depth); %kg
boat.added_mass_y = pi*rho*(boat.length/2)^2*(boat.depth); %kg
boat.Izz = 1; %kg*m
boat.added_Izz = rho*((boat.length/2)^2-(boat.width/2)^2)^2*boat.depth; %kg * m^2

%% "Aero" Params
switch aero_case
    case 'asym_point_particle'
        boat.cdx = 0.1;
        boat.cdy = 0.5;
        boat.k_drag_para = 1; 
        boat.k_drag_perp = 1; 
        boat.k_drag_rot = 0.4;
    case 'airfoil'
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

        boat.airfoil_lookup = [flip(aifoil_negative_side,1);airfoil_lookup];
        boat.alpha_lookup = airfoil_lookup(:,1);
        boat.cl_lookup = airfoil_lookup(:,2);
        boat.cd_lookup = airfoil_lookup(:,3);
end