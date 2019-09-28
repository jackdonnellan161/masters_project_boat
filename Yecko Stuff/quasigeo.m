% Barotropic Quasi-Geostrophic basin
% BT-QG model from Ocean Circulation Theory by Joe Pedlosky
% original code by Karl Helfrich, 2000; modified and extended by Phil Yecko, 2009;
%
% Pedlosky OCT Eqn (2.3.2):
% d(D^2(Psi))/dt + epsilon*J(Psi,D^2(Psi)) + dPsi/dx = wE - mu*D^2(Psi) + E*D^4(Psi)
% where D = del (nabla) operator
%      mu = ds   = R/(beta*L)
% epsilon = di^2 = U/(beta*L^2)
%       E = dm^3 = A_H/(beta*L^3)
% using scales:
%   L (horiz length), U (horiz velo), D (vert depth), Omega_0 (rotation)
% and with R=2*Omega_0*delta_E/(2D) 
% where delta_E is vert Ekman layer thickness

function out=quasigeo
    
 global X Y u v Psi dx dy
 global ib jb is js nx ny
 global mu epsilon E
 global Jacobian jpp jpx jxp
 global Px Py
 global time_now dt_plot mt_plot t_startsave t_starttracers
 global energy enstrophy t
 global fig_vel
 
 
% figure 1: u,v,Psi; special commands to make movies (comment in/out)
 fig_vel = figure(1);
 % set(fig_vel,'DoubleBuffer','on');
 % set(gca,'NextPlot','replace','Visible','off')
 % mov = avifile('floatie.avi')
 

% basin lengthscales 
Lx = 1.0;
Ly = 1.0;

% number of grid cells, ideally square and a power of 2 (for Poisson)
% ghost cells for no-slip BCs (Munk case)
nx = 128;
ny = 128;

% grid cell size (accounting for ghost zone)
% currently, must have dx = dy
dx = Lx / (nx-3);
dy = Ly / (ny-3);

% indicies of ghost, boundary, interior points:
ig = [1:nx];
jg = [1:ny];
ib = [2:nx-1];
jb = [2:ny-1];
is = [3:nx-2];
js = [3:ny-2];

% coordinate grids
x = dx * [-1:1:nx-2]';
y = dy * [-1:1:ny-2]';
[X,Y] = meshgrid (x,y);
X = X';
Y = Y';

% dt is timestep, mt is number of time steps
dt   = 0.01
tmin = 0.0;
tmax = 7.5;
mt   = round ((tmax-tmin)/dt);
it = [1:mt];
t = tmin + dt * [1:mt];

mt_plot = 10;   % how many steps between plots (1 unless dt<<1)
t_startsave=500 % need to adjust this ALWAYS based on model params
t_starttracers= mt+1 % after spin-up (500), use huge value for no tracers

% Stommel, Munk, and inertial params (see Pedlosky) normalized by L=Lx
ds = 0.04; mu=ds;
di = 0.02; epsilon=di^2; %0.02;
dm = 0.00; E=dm^3;
beta_eff = 1;

% boundary-condition flags (Munk case, E non-zero) 
% slip = 0 -> is no-slip (stress)
% slip = 1 -> is free-slip (no-stress);
slip_w = 0;
slip_e = 0;
slip_s = 0;
slip_n = 0;

% conserve enstrophy and energy
% set weights for Arakawa's J++, J+x and Jx+  
jpp =1/3;         % 1/3; 
jpx =1/3;         % 1/3; 
jxp =1-jpp-jpx;   % sum = 1 

% for periodic driving perturbations, amplitude and frequency
pertamp = 0.15;
omega   = 9.0;




% clean initialization
Psi       = zeros (nx, ny);
PsiStar   = zeros (nx, ny);
dPsidt    = zeros (nx, ny);
zeta_rel  = zeros (nx, ny);
u         = zeros (nx, ny);
v         = zeros (nx, ny);
f_RHS     = zeros (nx, ny);
Jacobian  = zeros (nx, ny);
energy    = zeros (mt,1);
enstrophy = zeros (mt,1);
w_Ekman   = zeros (nx,ny) ;

w_pert    = ones(nx,ny);

% initialize Ntrace initial tracer positions Px and Py
Ntrace      = 0;
Px          = zeros (Ntrace);
Py          = zeros (Ntrace);
%tracer_file = fopen('QG_out/talltest_3tracers.dat','w');




% initial ideal particle position(s) - need not fall on gridpoints:
% Consider: extend to inertial, Maxey-Riley <-- too slow, do post-run

% uninspired first try:
% Px(1) = 0.25;
% Py(1) = 0.25;

% Eric's 
%Px(1)=0.034
%Py(1)=0.266

% Eric's pair:
%Px(1)=0.142
%Py(1)=0.422
%Px(2)=0.142
%Py(2)=0.408

% 2nd pair: 26/7/09 run
%Px(1)=0.17
%Py(1)=0.394
%Px(2)=0.168
%Py(2)=0.446

% 3rd pair: 29/10/09 run
%Px(1)=0.162
%Py(1)=0.338
%Px(2)=0.162
%Py(2)=0.386

% 4th pair:
%Px(1)=0.161
%Py(1)=0.35
%Px(2)=0.162
%Py(2)=0.366

% 1st triple 12/3/2009:
%Px(1)=0.296
%Py(1)=0.248
%Px(2)=0.456
%Py(2)=0.248
%Px(3)=0.688
%Py(3)=0.248

   
point0  = zeros(2);
point1  = zeros(2);
point2  = zeros(2);
point3  = zeros(2);
s1      = zeros(2);
s2      = zeros(2);
s3      = zeros(2);
s4      = zeros(2);
delPxy  = zeros(2);





% time stepping 
% Huen method (RK2) for QG flow and RK4 for particles
% initial time, start stepping
time_now = tmin;

for tk = 1:1:mt
    
   time_now = t (tk);
   
   w_Ekman = zeros (nx,ny) ; % excess?  BCs?
   w_pert  = ones(nx,ny);    % excess?
   w_Ekman (is,js)=-sin(2*pi/Ly*Y(is,js)-pertamp*(2*pi/Ly)*sin(omega*time_now)*w_pert(is,js));


   % PREDICTOR

   % apply BCs
   % no normal flow:
   Psi (   2,:   ) = zeros (1 ,ny);   % western boundary
   Psi (nx-1,:   ) = zeros (1 ,ny);   % eastern "
   Psi (:   ,   2) = zeros (nx, 1);   % southern "
   Psi (:   ,ny-1) = zeros (nx, 1);   % northern "
   %  
   % tangential-flow (Munk case):
   % no-slip  : d(Psi)/dn  = 0,  ghost points equal to first interior point.
   % free-slip: del^2(Psi) = 0,  ghost points equal to negative of first interior point. 
   %
   Psi (1 ,: ) = (-1)^slip_w * Psi (   3,:   );   % western boundary
   Psi (nx,: ) = (-1)^slip_e * Psi (nx-2,:   );   % eastern "
   Psi (: , 1) = (-1)^slip_s * Psi (:   ,   3);   % southern " 
   Psi (: ,ny) = (-1)^slip_n * Psi (:   ,ny-2);   % northern "
   
   
   
   
   

   % RHS of elliptic eqn del^2(dPsidt) = f_RHS
   Del2Psi = zeros (nx,ny);
   Del4Psi = zeros (nx,ny);
   f_RHS   = zeros (nx,ny);

   % calculate del^2(Psi)
   % swap x and y for consistency with matlab sense of matrix orientation
   % factor 4 for consistency with matlab del2 (Laplacian)
   Del2Psi (ig,jg) = 4 * del2 (Psi (ig,jg), dy, dx);

   % relative vorticity is del^2 Psi
   zeta_rel (ig,jg) = Del2Psi;

   % calculate del^4 Psi for lateral friction term (Munk)
   if dm ~= 0
     Del4Psi (ig,jg) = 4*del2 (Del2Psi (ig,jg), dy, dx);
   end

   % compute velocity 
   velocity (Psi (ig,jg));

   % compute Jacobian
   comp_Jacobian (Psi (ig,jg), zeta_rel (ig,jg), u (ig,jg), v (ig,jg));

   % compute RHS forcing, mult by dx^2 (assumes dy=dx) in prep for elliptic solver.
   f_RHS (is,js)=dx^2*(w_Ekman (is,js)-beta_eff*v(is,js)-mu*Del2Psi(is,js) ...
                 -epsilon*Jacobian(is,js)+E*Del4Psi (is,js));

   
   % solve elliptic problem
   dPsidt (ib,jb) = poisson_solve (f_RHS (ib,jb)) ;
 
   % provisional PsiStar at new-time
   PsiStar (ib,jb) = Psi (ib,jb) +  dPsidt (ib,jb) ;

   
   
   % CORRECTOR: 
   
   
   % BCs ON PsiStar:
   % no normal-flow: 
   %
   PsiStar (   2,:   ) = zeros (1 ,ny);   % western boundary
   PsiStar (nx-1,:   ) = zeros (1 ,ny);   % eastern "
   PsiStar (:   ,   2) = zeros (nx, 1);   % southern "
   PsiStar (:   ,ny-1) = zeros (nx, 1);   % northern "
   %  
   % tangential-flow condition:
   % (only relevant for munk boundary layers)
   % no-slip : d(Psi)/dn  = 0,  ghost points equal to first interior point.
   % free-slip: del^2(Psi) = 0,  ghost points equal to negative of first interior point. 
   %
   PsiStar (1 ,: ) = (-1)^slip_w * PsiStar (   3,:   );   % western boundary
   PsiStar (nx,: ) = (-1)^slip_e * PsiStar (nx-2,:   );   % eastern "
   PsiStar (: , 1) = (-1)^slip_s * PsiStar (:   ,   3);   % southern " 
   PsiStar (: ,ny) = (-1)^slip_n * PsiStar (:   ,ny-2);   % northern "

   
   

   
   % evaluate RHS of elliptic eqn D^2(dPsidt) = f_RHS

   % initialize some arrays
   Del2Psi = zeros (nx,ny);
   Del4Psi = zeros (nx,ny);
   f_RHS = zeros (nx,ny);

   % calculate del^2 Psi
   % swap x and y for consistency with matlab sense of matrix orientation
   % factor 4 for consistency with matlab del2 (Laplacian)
   Del2Psi (ig,jg) = 4 * del2 (PsiStar (ig,jg), dy, dx);

   zeta_rel (ig,jg) = Del2Psi;

   % del^4 Psi for lateral friction term (Munk case)
   if dm ~= 0
      Del4Psi (ig,jg) = 4*del2 (Del2Psi (ig,jg), dy, dx);
   end

   velocity (PsiStar (ig,jg));

   % calculate Jacobian
   comp_Jacobian (PsiStar (ig,jg), zeta_rel (ig,jg), u (ig,jg), v (ig,jg));

   % compute RHS forcing, mult by dx^2 (assumes dy=dx) in prep for elliptic solver.
   f_RHS (is,js) = dx^2*(w_Ekman (is,js)-beta_eff*v (is,js)-mu*Del2Psi (is,js) ...
              - epsilon*Jacobian (is,js)+ E*Del4Psi (is,js));

   
   
   dPsidt (ib,jb) = poisson_solve (f_RHS (ib,jb));	      
   
   % new Psi:
   Psi (ib,jb) = 0.5*(Psi (ib,jb)+PsiStar(ib,jb))+0.5*dt*dPsidt(ib,jb);
           
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% update tracer(s) positions (aka compute pathlines):        
   
   if (tk>t_starttracers)
       
       % Loop over each tracer
       for kn=1:Ntrace

% easy but slow interpolation (Matlab built-in)
%        uPxinterp = griddata(X,Y,u,Px(kn),Py(kn));
%        vPyinterp = griddata(X,Y,v,Px(kn),Py(kn));


% fast, better bilinear interpolation, see notes PY 2010
        point0=[Px(kn) Py(kn)];
        s1=uvinterp(point0);
        point1=point0+(dt/2)*s1;
        s2=uvinterp(point1);
        point2=point0+(dt/2)*s2;
        s3=uvinterp(point2);
        point3=point0+dt*s3;
        s4=uvinterp(point3);
        
        % RK4 tracer step
        delPxy = (dt/6)*(s1+2*s2+2*s3+s4);

%        delP = norm(delPxy);        
%        if (delP > dx/2) 
%            fprintf(' TRACER STEP TOO BIG: %e VS %e \n',delP,dx);
%        end
                
        Px(kn) = Px(kn) + delPxy(1);
        Py(kn) = Py(kn) + delPxy(2);
       
       end % kn=1:Ntrace

%%% 3 particle:
       fprintf(tracer_file,'%e ',tk,Px(1),Py(1),Px(2),Py(2),Px(3),Py(3));

%% 2 particle       fprintf(tracer_file,'%e ',tk,Px(1),Py(1),Px(2),Py(2));
% 1 particle:       fprintf(tracer_file,'%e ',tk,Px(1),Py(1));
       
    fprintf(tracer_file,'\n');

   end % pathline section (tk>t_starttracers)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
           
           
   % compute kinetic energy and enstrophy (or not, comment out)
    energy (tk) = 0.5 * (sum(sum(u.^2))+sum(sum(v.^2))) / (nx*ny);
    enstrophy (tk) = 0.5 * sum(sum(zeta_rel.^2)) / (nx*ny);             
   
  % if (tk>t_starttracers)
     if (mod (tk,mt_plot) == 0) 

        plot_velocity;

        %%% MOVIE FRAMES        
        framefile = ['QG_movie/frame_' num2str(tk,'%05d') '.png' ]; 
        print ('-dpng',framefile,'-r0');     

       %% save results to disk file named QGrun_no.mat
       %file_save = ['QG_runs/QG_run_' num2str(run_no) '_time_' num2str(tk) '.mat' ];
       %save (file_save,'time_now', 'Psi', 'zeta_rel', 'u', 'v');
       
       
     end     
  % end
 
   
   
% save data after spin-up time (t>t_startsave):
% save results to disk file named <<SOMETHING>>.mat
   if (tk>t_startsave) 
%    file_data = ['QG_out/pert0p1_time_' num2str(tk) '.mat' ];
%    file_data = ['QG_out/pa15om1ds04di02_g128t01_time_' num2str(tk,'%05d') '.mat' ];
     file_data = ['QG_out/pa15om9stom_04_02_g128t001_time_' num2str(tk,'%05d') '.mat' ];

       % also Psi save option
       %%file_data = ['QG_out/Psi_time_' num2str(tk,'%05d') '.mat' ];
       %%save ( file_data, 'Psi');     
     save ( file_data, 'u', 'v');     
   end
  
   if (tk==1) 
       % save results to disk file named onelay_time_xxx.mat
       file_grid = ['QG_out/grid128.mat' ];
       save ( file_grid, 'X', 'Y');
   end
   

end   % tk=1:mt


% figure 2: KE(t),Ens(t); if not wanted: comment out
 plot_en_ens;

return
end  % end of main function: quasigeo 
   
   

% AUXILIARY FUNCTIONS: uvinterp, comp_Jacobian, velocity,
% poisson_solve, plot_velocity, plot_en_ens


% Interpolation of flow velocity for (ideal) particles

function uivi = uvinterp(p)
global X Y u v dx dy

veli = floor(p(1)/dx)+2;
velj = floor(p(2)/dy)+2;
        
xrel = p(1)-X(veli,velj);
yrel = p(2)-Y(veli,velj);
       
% BILINEAR interpolation routine, see PY notes 2010
   uPxinterp = u(veli,    velj)  * (dx-xrel)*(dy-yrel) ...
             + u(veli+1,  velj)  * xrel*(dy-yrel) ...
             + u(veli,  velj+1)  * (dx-xrel)*yrel ...
             + u(veli+1,velj+1)  * xrel*yrel;
   uivi(1) = uPxinterp / (dx*dy);
        
   vPyinterp = v(veli,velj)*(dx-xrel)*(dy-yrel) ...
             + v(veli+1,velj)*xrel*(dy-yrel) ...
             + v(veli,velj+1)*(dx-xrel)*yrel ...
             + v(veli+1,velj+1)*xrel*yrel;
   uivi(2) = vPyinterp/(dx*dy);
   % end of function uvinterp
   
   return
   end

   %%%%%%%%%
   
   function comp_Jacobian (Psi_in, zeta_in, u_in, v_in) ;

% Compute Arakawa Jacobain term of non-linear advection
% see Arakawa (1966,1997) or Haltiner & Williams (1980)

   % WARNING: extended BUT UNTESTED to dy NOT EQUAL dx

   % an alternate Jacobian to consider: Sadourny (1975)
    
global is js nx ny dx dy
global Jacobian jpp jpx jxp

jacob_pp = zeros (nx,ny);
jacob_xp = zeros (nx,ny);
jacob_px = zeros (nx,ny);
                      
%  J++ AKA -J1 
if jpp~=0
    jacob_pp (is, js) =((Psi_in  (is+1,js  ) - Psi_in  (is-1,js  ))   ...
                     .* (zeta_in (is  ,js+1) - zeta_in (is  ,js-1))   ...
                      - (Psi_in  (is  ,js+1) - Psi_in  (is  ,js-1))   ...
                     .* (zeta_in (is+1,js  ) - zeta_in (is-1,js  )))  ...
                      / (4*dx*dy);
end
        
%  J+x
if jpx~=0
    jacob_px (is, js) =(Psi_in  (is+1,js  )                          ...
                     .*(zeta_in (is+1,js+1) - zeta_in (is+1,js-1))   ...
                      - Psi_in  (is-1,js  )                          ...
                     .*(zeta_in (is-1,js+1) - zeta_in (is-1,js-1))   ...
                      - Psi_in  (is  ,js+1)                          ...
                     .*(zeta_in (is+1,js+1) - zeta_in (is-1,js+1))   ...
                      + Psi_in  (is  ,js-1)                          ...
                     .*(zeta_in (is+1,js-1) - zeta_in (is-1,js-1)))  ...
                     / (4*dx*dy);
end
            
%  Jx+
if jxp~=0
   jacob_xp (is,js) =(zeta_in (is  ,js+1)                     ...
                   .*(Psi_in  (is+1,js+1) - Psi_in (is-1,js+1))  ...
                    - zeta_in (is  ,js-1)                     ...
                   .*(Psi_in  (is+1,js-1) - Psi_in (is-1,js-1))  ...
                    - zeta_in (is+1,js  )                     ...
                   .*(Psi_in  (is+1,js+1) - Psi_in (is+1,js-1))  ...
                    + zeta_in (is-1,js  )                     ...
                   .*(Psi_in (is-1,js+1) - Psi_in (is-1,js-1))) ...
                    /(4*dx^2);
end

% the famous conserving combo:
Jacobian (is,js) = jpp * jacob_pp (is,js) ...
                 + jpx * jacob_px (is,js) ...
                 + jxp * jacob_xp (is,js) ; 

return
end

%%%%%%%


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

%%%%%%%

function u = poisson_solve (p)

% modified fast Poisson solver with boundary values, original by Chick Denham
% for more info & updates see:
% WHOI Matlab Snack Bar: woodshole.er.usgs.gov/staffpages/cdenham/public_html/snackbar/snackbar.html
% or GitHub sea-mat: github.com/sea-mat/seagrid/blob/master/seagrid/fps.m

q = p;

% extract the interior.
q = q (2:end-1, 2:end-1);

% odd-symmetry, sine-transform scheme
theSign = -1;   

[m, n] = size(q);
q = [zeros(m, 1) q zeros(m, 1) theSign*fliplr(q)];
q = [zeros(1, 2*n+2); q; zeros(1, 2*n+2); theSign*flipud(q)];

% fast poisson transform on square cell grid 
% (dx^2 factor already multiplies p).

[mm, nn] = size (q);
i = (0:mm-1).' * ones(1, nn);
j = ones (mm, 1) * (0:nn-1);
weights = 2 * (cos (2*pi*i/mm) + cos (2*pi*j/nn) - 2);
weights (1, 1) = 1;
pp = ifft2 (fft2 (q) ./ weights);
if isreal (q), pp = real(pp); end

% retain relevant piece.
u = p;
u (2:end-1, 2:end-1) = pp (2:m+1, 2:n+1);

return
end



%%%%%%%



function plot_velocity

global X Y u v Psi
global ib jb
global Px Py time_now
global fig_vel fig_en_ens

% axes ignoring ghost layers
Lxi = X (2    ,2    );
Lxf = X (end-1,2    );
Lyi = Y (2    ,2    );
Lyf = Y (2    ,end-1);

KE = sqrt (u(:,:).*u(:,:) + v(:,:).*v(:,:));
KE_max = max(max(KE));

figure (fig_vel);
clf;

% qstag (integer) often needed to re-space the (u,v) vectors!
qstag = 4;
qib=qstag*(1:length(ib)/qstag);
qjb=qstag*(1:length(jb)/qstag);

% plot surface
quiver (X (qib,qjb), Y (qib,qjb), u (qib,qjb), v (qib,qjb));
hold on;
contour (X (ib,jb), Y (ib,jb), Psi (ib,jb));
%colorbar; % if desired, else comment out

% Need to generalize below to N tracers!
%plot(Px(1),Py(1),'k*','MarkerSize',10);
%plot(Px(2),Py(2),'ro','MarkerSize',10);
%plot(Px(3),Py(3),'gd','MarkerSize',10);

axis ([Lxi Lxf Lyi Lyf]);
axis equal ;

% annotation
xlabel('x');
ylabel('y');
%zlabel('velocity')
%title ([' Time = ' num2str(time_now,2) '   Max(velo) = ' num2str(KE_max,2)]);
title(['velocity and streamfunction at time = ' num2str(time_now)]);

drawnow;

return
end

%%%%%%%%%%%

function plot_en_ens

% plot domain-integrated kinetic energy and enstrophy timeseries

global energy enstrophy
global t mt time_now

figure(2);

% energy
% subplot (2,1,1);
plot (t,energy);

title (['total energy vs time']);
xlabel ('time '); 
ylabel ('energy ');

% enstrophy
% subplot (2,1,2);
% plot (t,enstrophy)
% xlabel ('time '); 
% ylabel ('enstrophy');

%drawnow;

return
end
