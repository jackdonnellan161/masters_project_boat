function [U_field,V_field,U_field_load,V_field_load] = velocity_field_init(field_str, param_struct)
    U_field = zeros(length(param_struct.xs),length(param_struct.ys));
    V_field = zeros(length(param_struct.xs),length(param_struct.ys));
    U_field_load = zeros(length(param_struct.xs),length(param_struct.ys));
    V_field_load = zeros(length(param_struct.xs),length(param_struct.ys));
    switch field_str
        case 'no_flow'
            
        case 'gyre'
            f = param_struct.a.*param_struct.xs.^2 + param_struct.b.*param_struct.xs;
            dfdx = (2.*param_struct.a.*param_struct.xs + param_struct.b)';
            U_field = (-pi.*param_struct.A_field.*sin(pi.*f)'*cos(pi.*param_struct.ys))'; %m/s
            V_field = (pi.*param_struct.A_field.*dfdx.*cos(pi.*f)'*sin(pi.*param_struct.ys))'; %m/s
            
        case 'QG_flow'
            
            
            
            
    end
            
    U_field_load = flipud(U_field);
    V_field_load = flipud(V_field);
end