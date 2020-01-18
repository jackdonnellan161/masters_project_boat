%% John Donnellan
%  Consolidated boat model parameters with switch for different aero models

function boat = boat_init(aero_case,rho)
%% Size/Dimensions
boat.length = 0.128; %m %From Michini Thesis
boat.width = 0.063; %m %From Michini Thesis
boat.depth = 0.032; %m %From Michini Thesis
boat.d_motors = boat.width-0.02; %distance between the centers of motor shafts, m
boat.Ax = boat.width*boat.depth; %m^2
boat.Ay = boat.depth*boat.length; %m^2

%% Mass/Inertias
boat.m = 0.1095; %kg %From Michini Thesis
boat.added_mass_x = pi*rho*(boat.width/2)^2*(boat.depth); %kg
boat.added_mass_y = boat.added_mass_x;
% boat.added_mass_y = pi*rho*(boat.length/2)^2*(boat.depth); %kg
boat.Izz = 1.2e-4; %kg*m^2 %From Michini Thesis
boat.added_Izz = rho*((boat.length/2)^2-(boat.width/2)^2)^2*boat.depth; %kg * m^2

%% "Aero" Params
switch aero_case
    case 'asym_point_particle'
        boat.cdx = 1;
        boat.cdy = 1;
        boat.k_drag_para = 8*boat.width; %From Michini Thesis
        boat.k_drag_perp = 1.28; %Flat Plate Drag
        boat.k_drag_rot = 1; %From Michini Thesis
    case 'airfoil'
        boat.airfoil_data = csvread('airfoil_data.csv');
        boat.Alpha_lookup = boat.airfoil_data(:,1);
        boat.CL_lookup = boat.airfoil_data(:,2);
        boat.CD_lookup = boat.airfoil_data(:,3);
        boat.CN_lookup = boat.airfoil_data(:,4);
        boat.CA_lookup = boat.airfoil_data(:,5);
end

%% Prop System
boat.Vmax = 0.2; %m/s
%boat.Tmax = 0.5*(0.5*rho*boat.Vmax^2*boat.k_drag_para*boat.Ax); %N 
%Tmax for non-linear drag model, current model uses linear drag so this
%thrust is insufficient by ~ a factor of 4
boat.Tmax = boat.Vmax.*boat.k_drag_para/2;