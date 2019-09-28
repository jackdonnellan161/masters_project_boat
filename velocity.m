function velocity (Psi_in)

% compute flow velocities from streamfun: u = -dPsi/dy, v = dPsi/dx

global u v
global ib jb nx ny dx dy

% initialize to zero -- needed?
u = zeros (nx,ny);
v = zeros (nx,ny);

u (: ,jb) = -(Psi_in (:,jb+1) - Psi_in (:,jb-1)) / (2*dy);
v (ib,: ) =  (Psi_in (ib+1,:) - Psi_in (ib-1,:)) / (2*dx);

return
end
