%% Discrete Elastic Rods (Bergou et al, SIGGRAPH, 2010) implementation
% Khalid Jawed, khalidjm@mit.edu
% May 2016
% A very fast C++ implementation can be found at
% http://www.cs.columbia.edu/cg/elastic_coiling/
% This MATLAB code is for demonstration only, and is very slow

%%
clear all;
clc;

fprintf('Discrete elastic rods\n');

%%
global m EI EA GJ dt x0 u nv ne d1 d2 m1 m2 x refLen voronoiRefLen tangent
global xCons nCons uUncons garr ctime RodLength visc r0
global consInd unconsInd
global ScaleSolver theta0 refTwist kappaBar
global tol maximum_iter
global Y G nu
%% Inputs
% number of vertices
nv = 50;

% Time step
dt = 1e-2;

% Rod Length
RodLength = 0.1;

% Natural radius of curvature of the rod
R = 0;

% Density
rho = 1000;

% Cross-sectional radius of rod
r0 = 1e-3;

% Young's modulus
Y = 1e6;

% Poisson ratio
nu = 0.5;

% Shear modulus
G = Y/(2.0*(1.0+nu));

% gravity
% g = [0, 0, -9.81];
g = [0, 0, 0];

% viscosity
visc = 0; % doesn't matter

% Tolerance on force function. This is multiplied by ScaleSolver so that we
% do not have to update it based on edge length and time step size
tol = 1e-7;

% Maximum number of iterations in Newton Solver
maximum_iter = 100;

% Total simulation time (it exits after t=totalTime)
% totalTime = 50;
totalTime = 50;

% Indicate whether images should be saved
saveImage = 1;

% How often the plot should be saved? (Set plotStep to 1 for each plot to
% be saved)
% plotStep = 10;
plotStep = 10;
%% Utility quantities
ne = nv - 1;
EI = Y * pi * r0^4 / 4;
EA = Y * pi * r0^2;
GJ = G * pi * r0^4/2;
dm = pi * r0^2 * RodLength * rho / ne;

%% Geometry of the rod
radius = 0.1*0.1/(2*pi); %helix radius;
pitch = 0.01; %helix pitch;
T_prep = 10; %preparation time, before forming a helix;

b_val = pitch/(2*pi); %b-constant;
t_val = RodLength/sqrt(radius^2 + b_val^2);
torsion = b_val/(radius^2 + b_val^2);
endTheta = -torsion * RodLength;
slope = b_val/radius; %helix slope;
del_l = RodLength/ne;
del_t = t_val/ne; 

%% Initial Orientation as Helix
%Helical shaped rod: position of the first two nodes
first1X = 0; first2X = b_val * del_t;
first1Y = 0; first2Y = radius * sin(del_t);
first1Z = radius; first2Z = radius * cos(del_t);
%Helical shaped rod: position of the last two nodes
last1X = b_val * ne * del_t; last2X = b_val * (ne-1) * del_t;
last1Y = radius * sin(ne * del_t); last2Y = radius * sin((ne-1) * del_t);
last1Z = radius * cos(ne * del_t); last2Z = radius * cos((ne-1) * del_t);
%Helical shaped rod: twist


%% Geometry of the rod
nodes = zeros(nv, 3);
if (R==0)
    for c=1:nv
%         tt = (c-1) * del_t;
%         nodes(c, 1) = b_val * tt; %%x-coord,x = bt;
%         nodes(c, 2) = radius * sin(tt); %%y-coord,y = asin(t);
%         nodes(c, 3) = radius * cos(tt); %%z-coord,z = acos(t);
           nodes(c, 1) = (c-1) * RodLength / ne;
    end
else
    endTheta = RodLength / R;
    for c=1:nv
        thTmp = (c-1) * endTheta / ne;
        nodes(c, 1) = R - R * cos(thTmp);
        nodes(c, 2) = R * sin((c-1) * endTheta / ne);
    end
end

%% Multiplier for Force & Jacobian
ScaleSolver = dm *(RodLength/ne) /dt^2;

%% Compute Mass
m = zeros(3*nv+ne, 1);
for c=1:nv
    if c==1
        m( 4 * (c-1) + 1 : 4 * (c-1) + 3) = dm/2;
    elseif c==nv
        m( 4 * (c-1) + 1 : 4 * (c-1) + 3) = dm/2;
    else
        m( 4 * (c-1) + 1 : 4 * (c-1) + 3) = dm;
    end
end
for c=1:ne
    m( 4 * c ) = dm/2 * r0^2; % I = 1/2 m r ^ 2
end

%% gravity
garr = zeros(3*nv+ne, 1);
for c=1:nv
    garr( 4 * (c-1) + 1 : 4 * (c-1) + 3) = g;
end

%% Reference length and Voronoi length
refLen = zeros(ne, 1);
for c=1:ne
    dx = nodes(c+1, :) - nodes(c, :);
    refLen(c) = norm(dx);
end
voronoiRefLen = zeros(nv, 1);
for c=1:nv
    if c==1
        voronoiRefLen(c) = 0.5 * refLen(c);
    elseif c==nv
        voronoiRefLen(c) = 0.5 * refLen(c-1);
    else
        voronoiRefLen(c) = 0.5 * (refLen(c-1) + refLen(c));
    end
end

%% Reference director & material director
d1 = zeros(ne, 3); % reference director, u (or d1)
d2 = zeros(ne, 3); % reference director, v (or d2)
tangent = zeros(ne, 3); % tangent
for c=1:ne
    dx = nodes(c+1,:) - nodes(c,:);
    tangent(c,:) = dx / norm(dx);
end

% Figure out a good choice for d1(1)
t0 = tangent(1,:);
t1 = [0 0 -1];
d1Tmp = cross(t0, t1);
if (abs(d1Tmp) < 1.0e-6)
    t1 = [0 1 0];
    d1Tmp = cross(t0, t1);
end
d1(1,:) = d1Tmp;

d1_l = d1(1,:);
d2(1,:) = cross(t0, d1_l);
for c=2:ne
    t1 = tangent(c,:);
    d1_l = parallel_transport(d1_l, t0, t1);
    d1_l = (d1_l - dot(d1_l, t1) * t1);
    d1_l = d1_l / norm(d1_l);
    d1(c,:) = d1_l;
    d2_l = cross(t1, d1_l);
    d2(c,:) = d2_l;
    t0 = t1;
end

%% Initial
x0 = zeros(3*nv + ne, 1);
for c=1:nv
    x0( 4 * (c-1) + 1) = nodes(c,1);
    x0( 4 * (c-1) + 2) = nodes(c,2);
    x0( 4 * (c-1) + 3) = nodes(c,3);
end
x0(4:4:end) = 0; % theta
x = x0;

%% Constrained dofs
dummyInd = 1:length(x);
consInd = [1:7, 3*nv + ne-6: 3*nv + ne]; % constrained dof
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

xCons = x(consInd); % first two nodes and one edge angle
nCons = length(xCons);
u = (x - x0) / dt;
uUncons = u(unconsInd); % unconstrained dof

%% Compute material director
theta0 = zeros(ne,1);
[m1, m2] = computeMaterialDirectors(d1, d2, theta0);
refTwist = zeros(ne,1);
% refTwist = getRefTwist(d1, tangent, refTwist); % Should be all zero

%% Natural curvature computation
kappaBar = getkappaBar( x, m1, m2 );

%% Create director to save image
imageDirectory = date;
if (saveImage~=0)
    mkdir(imageDirectory);
end

%% Options for MATLAB optimization
% This is only used if simple Newton's method fails.
% We increase MaxFunEvals to 10 * maximum_iter
options = optimoptions(@fsolve, 'Display','iter', ...
    'Algorithm', 'trust-region-reflective', ...
    'Jacobian','off', ...
    'TolFun', tol, 'TolX', tol, ...
    'MaxFunEvals', maximum_iter*10, 'MaxIter', maximum_iter*10);

%% Time marching
Nsteps = round(totalTime/dt); % number of time steps
Tsteps = round(T_prep/dt); %number of time steps before it forms a helix;

%% Twist Analysis
allAngle = zeros(Nsteps, ne);
allTwistAngle = zeros(Nsteps, ne);
x4angle = zeros(Tsteps,1);

ctime = 0;

rotateAngle = zeros(Nsteps,1);
normU = zeros(Nsteps,1);
error_from_helix = zeros(Nsteps,1);
%  finishDrag = 0;
R0 = 0; P0 = 1000*RodLength; %make it large;
for timeStep=1: Nsteps
    
    fprintf('t=%f\n', ctime);
    
    xUncons = x(unconsInd) + uUncons * dt; % Guess
    
if (timeStep <= Tsteps)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Geometry of the rod
tempR = R0 + (timeStep/Tsteps)*(radius - R0); %helix radius;
tempP = P0 + (timeStep/Tsteps)*(pitch - P0); %helix pitch;

tempB = tempP/(2*pi); %b-constant;
tempT = RodLength/sqrt(tempR^2 + tempB^2);
tempTor = tempB/(tempR^2 + tempB^2);
tempendTheta = -tempTor * RodLength;
temp_del_l = RodLength/ne;
temp_del_t = tempT/ne; 

%% Initial Orientation as Helix
%Helical shaped rod: position of the first two nodes
first11X = 0; first22X = tempB * temp_del_t;
first11Y = 0; first22Y = tempR * sin(temp_del_t);
first11Z = tempR; first22Z = tempR * cos(temp_del_t);
%Helical shaped rod: position of the last two nodes
last11X = tempB * ne * temp_del_t; last22X = tempB * (ne-1) * temp_del_t;
last11Y = tempR * sin(ne * temp_del_t); last22Y = tempR * sin((ne-1) * temp_del_t);
last11Z = tempR * cos(ne * temp_del_t); last22Z = tempR * cos((ne-1) * temp_del_t);

%Last node;
        u1x = (x(4*nv-3)-last11X); 
        u1y = (x(4*nv-2)-last11Y); 
        u1z = (x(4*nv-1)-last11Z); 
        v1 = [u1x, u1y, u1z]'; 
        x(3*nv + ne-2:3*nv + ne) = x(3*nv + ne-2:3*nv + ne) - v1 * 1 ;
%Second to the Last node;
        u2x = (x(4*nv-7)-last22X);
        u2y = (x(4*nv-6)-last22Y);
        u2z = (x(4*nv-5)-last22Z);
        v2 = [u2x, u2y, u2z]';
        x(3*nv + ne-6:3*nv + ne-4) = x(3*nv + ne-6:3*nv + ne-4) - v2 * 1 ;
%First node;
        u11x = (x(1)-first11X); 
        u11y = (x(2)-first11Y); 
        u11z = (x(3)-first11Z); 
        v11 = [u11x, u11y, u11z]'; 
        x(1:3) = x(1:3) - v11 * 1 ;
%Second node;
        u22x = (x(5)-first22X);
        u22y = (x(6)-first22Y);
        u22z = (x(7)-first22Z);
        v22 = [u22x, u22y, u22z]';
        x(5:7) = x(5:7) - v22 * 1 ;
%Temp angle;
        x(4) = tempendTheta + sum(refTwist); 
%          x(4) = 0 + (timeStep/Tsteps)*(endTheta - 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     flag = 0;
%     tolX = (RodLength / ne)/5;
%     if (sqrt((x(4*nv-3)-last1X)^2+(x(4*nv-2)-last1Y)^2+(x(4*nv-1)-last1Z)^2)>=tolX)
%         u1x = (x(4*nv-3)-last1X); 
%         u1y = (x(4*nv-2)-last1Y); 
%         u1z = (x(4*nv-1)-last1Z); 
%         v1 = [u1x, u1y, u1z]'; 
%         x(3*nv + ne-2:3*nv + ne) = x(3*nv + ne-2:3*nv + ne) - v1 * dt;
% %         flag = 1;
%     end
%     if(sqrt((x(4*nv-7)-last2X)^2+(x(4*nv-6)-last2Y)^2+(x(4*nv-5)-last2Z)^2)>=tolX)
%         u2x = (x(4*nv-7)-last2X);
%         u2y = (x(4*nv-6)-last2Y);
%         u2z = (x(4*nv-5)-last2Z);
%         v2 = [u2x, u2y, u2z]';
%         x(3*nv + ne-6:3*nv + ne-4) = x(3*nv + ne-6:3*nv + ne-4) - v2 * dt;
%         flag = 1;
%     end
%     if(sqrt((x(1)-first1X)^2+(x(2)-first1Y)^2+(x(3)-first1Z)^2)>=tolX)
%         u11x = (x(1)-first1X); 
%         u11y = (x(2)-first1Y); 
%         u11z = (x(3)-first1Z); 
%         v11 = [u11x, u11y, u11z]'; 
%         x(1:3) = x(1:3) - v11 * dt;
%         flag = 1;
%     end
%     if(sqrt((x(5)-first2X)^2+(x(6)-first2Y)^2+(x(7)-first2Z)^2)>=tolX)
%         u22x = (x(5)-first2X);
%         u22y = (x(6)-first2Y);
%         u22z = (x(7)-first2Z);
%         v22 = [u22x, u22y, u22z]';
%         x(5:7) = x(5:7) - v22 * dt;
%         flag = 1;
%     end
%     if flag == 0
%         finishDrag = 1;
%     end
else %start twisting;
        omega = 3*pi/10;
        x(4) = x(4) -  omega * dt;
end
    
    xCons = x(consInd);
%  % Calculate the difference 
%     rotateAngle(timeStep) = x(4);
    
    % Attempt 1 for solution
    [xUncons, error] = objfun(xUncons);
    
    x(consInd) = xCons;
    x(unconsInd) = xUncons;
    u = (x - x0) / dt; % velocity
    ctime = ctime + dt; % current time
    uUncons = u(unconsInd);
    
    % Theta
    theta0 = x(4:4:end);
    
    % Update material and reference directors
    tangent = computeTangent(x);
    [d1, d2] = computeTimeParallel(d1, x0, x);
    refTwist = getRefTwist(d1, tangent, refTwist); % Compute reference twist
    [m1, m2] = computeMaterialDirectors(d1, d2, theta0);
    
%     totdiff = zeros(nv,1);
%     for count = 1 : nv
%         dxx = x(4 * (count-1) + 1);
%         dyy = x(4 * (count-1) + 2);
%         dzz = x(4 * (count-1) + 3);
%         totdiff(count) = sqrt(dxx^2 + dyy^2 + dzz^2);
%     end
%     disDiff(timeStep) = mean(totdiff);
    % Update x0
%     normU(timeStep) = norm(u);
     % Calculate the difference 
    rotateAngle(timeStep) = x(4);
    error_abs = 0;
    for i = 1:nv
        temp_e = abs(sqrt(x(4*i-2)^2+x(4*i-1)^2) - tempR);
        error_abs = error_abs + temp_e;
    end
    error_avg = error_abs/nv;
    error_from_helix(timeStep) = error_avg;
    x0 = x;
    
    %% angle
%     for 
%     

for c = 1:ne
    if (c == 1)
        allTwistAngle(timeStep,c) = x(4*(c)) - 0 + refTwist(c);
    else
        allTwistAngle(timeStep,c) = x(4*(c)) - x(4*(c-1)) + refTwist(c);
    end
end

for c = 1:ne
    allAngle(timeStep,c) = x(4*c);
end
    
x4angle(timeStep,1) = x(4);
    %%
    
    if (mod(timeStep, plotStep) ==0)
        plotrod(x, d1, d2, m1, m2, ctime, saveImage, imageDirectory);
    end
end

% figure (2)
% plot(rotateAngle,normU);

re_angle = rotateAngle(Tsteps:Nsteps)-endTheta;
re_error = abs(error_from_helix(Tsteps:Nsteps));
% re_angle = rotateAngle(1:Tsteps);
% re_error = abs(error_from_helix(1:Tsteps));
figure (3)
plot(re_angle,re_error);
xlabel('Angle of rotation (rad)') 
ylabel('Deviation from Helix (m)') 
title('Helix Solution: Stability Test: Case 1. Negative Torsion')
grid on;

figure (4)
plot (x4angle);