clear all; clc;
close all; 

%% # grid particles and integration step/time
nx = 2^9; ny = nx/2;
dt = 0.1; t0 = 0; te = 15;
% double gyre params
params.epsilon = 1; params.A = 1; params.omega = pi;
lambda_bound = 0.3;

%% Set up grid with particles (in 2 dimensions)
x = linspace(0,2,nx); y = linspace(0,1,ny);
[xx,yy] = meshgrid(x,y);
X0 = [xx(:), yy(:)];

%% Integrate particles for some time T
tic
h = dt; T = te - t0; n = T/h;
Xe = rk4t(@(t,x) doublegyre(t,x,params),X0,h,n,t0);

%% Compute FTLE field
toc
lambda = ftle(X0,Xe,nx,ny,T);
toc

%% Plot FTLE field
close all
figure, pcolor(xx,yy,lambda), shading flat, colorbar
caxis([0 max(max(lambda))])
title(['Forward FTLE field of double gyre, T = ',num2str(T)] )
axis('equal')

map = [1, 0, 0 
    1, 1, 1
    0, 0, 1];

lambda(lambda<lambda_bound) = 0; 
figure(2), hold on
pcolor(xx,yy,lambda)
colormap(map), caxis([-max(max(lambda)) max(max(lambda))])
shading flat, colorbar
title(['Forward FTLE field of double gyre: ',...
      '[', num2str(t0), ', ', num2str(te),']'])
axis('equal')

lambdaf = lambda;