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
        boat.airfoil_data = csvread('airfoil_data.csv');
        boat.Alpha_lookup = boat.airfoil_data(:,1);
        boat.CL_lookup = boat.airfoil_data(:,2);
        boat.CD_lookup = boat.airfoil_data(:,3);
        boat.CN_lookup = boat.airfoil_data(:,4);
        boat.CA_lookup = boat.airfoil_data(:,5);
end