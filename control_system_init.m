%% John Donnellan
%  Consolidated control system parameters with switch for different models

function [controls, actuators, nav] = control_system_init(controls_case, actuators_case, nav_case)
    %% Control Params
    switch controls_case
        case 'PI'
            controls.kpx = 1;
            controls.kpy = 1;
            controls.kix = 0.1;
            controls.kiy = 0.1;
        case 'LQR'
    end

    %% Actuator Params
    actuators.L_motor = 0.1;
    actuators.kphi_motor = 0.3;
    actuators.J_motor = 0.1;
    actuators.b_motor = 0.01;
    actuators.R_motor = 2.0;
    switch actuators_case
        case 'use_motors'
            actuators.enable = 1;
        case 'direct_thrust'
            actuators.enable = 0;
    end
    
    %% Nav Params
    switch nav_case
        case 'use_sensors'
            nav.sensors_enable = 1;
        case 'nav_passthrough'
            nav.sensors_enable = 0;
    end
end
