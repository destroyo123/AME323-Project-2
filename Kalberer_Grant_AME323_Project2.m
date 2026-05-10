%% AME 323 - Gas Dynamics, Project 2: Wind-Tunnel Design
% *Authors:* Maren Kalberer, Etan Grant
%%
% *Due Date:* 5-13-26
%% Description:
% 
% Todo: write this later I guess?
% 
% 
% 
% 
% 

%% Housekeeping:
% Requires having the a few toolboxes (Aerospace Toolbox) installed for certain functions

% Clear old variable values, outputs, and close old figures
close all;
clear;
clc;

% UofA colors!
arizonaBlue = [12, 35, 75] ./ 255;
arizonaRed = [171, 5, 32] ./ 255;
midnight = [0, 28, 72] ./ 255;
azurite = [30,82,136] ./ 255;
oasis = [55, 141, 189] ./ 255;
chili = [139, 0, 21] ./ 255;
bloom = [239, 64, 86] ./255;
sky = [192,211,235] ./255;
leaf = [112, 184, 101] ./255;
river = [0, 125, 132] ./255;
mesa = [169, 92, 66] ./255;

% Set interpereter to Latex. Lets us use subscripts and greek characters.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Given Values:
% Known relationships, and input parameters.

% Fluid properties
R_u = 8.314; % Universal Gas Const

% Molar masses
M_air = 28.965 / 1000; % kg/mol
M_Helium = 4.0026 / 1000; % kg/mol

M_fluid = M_air; % Since we're using air.

% Specific Gas Constants 
R = R_u/M_fluid; % (we're using air)

% Specific heat ratios of air and helium
gamma_air = 1.4;
gamma_helium = 1.67;

gamma = gamma_air; % since we're using air

% Liquid temperatures
Tl_helium = 2.2; % Kelvin
Tl_air = 50; % Kelvin

% Temperature the fluid MUST stay above
Tl = Tl_air; % Temperature our fluid (air) liquifies at

% Parameters

M_e = 6; % Exit mach
M_t = 1; % Throat

N = 4; %Number of waves

% h_n is given as the full diameter...
h_n = 350/2; %Nozzle height in mm
h_n = h_n/1000; % Converted to m


%% Functions
% Equations and Relations I mostly made for other homeworks & discussion
% problems. Equations come from NACA1135 unless otherwise stated.

% Area ratio equation
AeAt = @(M) (1/M) * (2/(gamma+1) * (1 + (((gamma-1)/2) * M^2 )))^( (gamma+1)/(2*(gamma-1)) );

%%
% *Prandtl-Meyer angle 'nu'*
% 
% *Inputs:*
%
%
% * Unperturbed Mach (AKA M1)
% 
% *Outputs:*
%
% * Prandtl-Meyer angle 'nu' (degrees)
%
    nu = @(M) sqrt((gamma+1)/(gamma-1)) * atand(sqrt(((gamma-1)/(gamma+1))*(M.^2 - 1))) - atand(sqrt(M.^2 - 1));
    % OUTPUTS DEGREES
%
%%
% *Mach from prandtl-meyer angle (IN DEGREES)* - Requires the Aerospace
% Toolbox
%
% *Inputs:*
%
% * Prandtl-Meyer angle of local flow 'nu' (degrees) (CAN BE AN ARRAY)
% * Flow ratio of specific heats 'gamma,' usually 1.4
%
% 
%
% *Outputs:*
%
% * Local mach number 'M'
%
%
function M = meyerMach(nu,gamma)
% * USES THE AEROSPACE TOOLBOX 'flowprandtlmeyer()
if(isscalar(nu))
    % If nu is a scalar, you can just use this function.
    M = flowprandtlmeyer(gamma, nu, 'nu');
else
    % if nu is a vector/matrix, you have to go 1 by 1 with this function
    M = zeros(length(nu), 1);
    for i = 1:length(nu)
        M(i,1) = flowprandtlmeyer(gamma, nu(i), 'nu');
    end
end
end
%%
% Post-Intersection nu
% 
% *Inputs:*
%
%
% * Left characteristic L = nu - delta
% * Right characteristic R = nu + delta
% 
% *Outputs:*
%
% * Post-Wave intersection Prandtl-Meyer angle 'nu' (degrees)
%


%% Pre-Calculations (Task 1)
% Easily-Derivable values found before calculating the complicated contour

% Ae/At also called (Ae/A*)
AreaRatio = AeAt(M_e);

% % Height of throat for a 3D nozzle with a circular cross-section.
% NOT USED FOR THIS PROJECT.
% h_t = (h_n ^2 /AreaRatio )^0.5;

% It was specified to do a planar / 2D nozzle... Like so:
% Height of throat for a 2D nozzle (rectangular cross section)
h_t = h_n/AreaRatio;

% Prandtl-Meyer angle 'nu' at exit
nu_e = nu(M_e); % degrees
nu_t = nu(M_t); % Degrees

% Max flow angle anywhere in the nozzle
delta_max = nu_e / 2; % degrees

% Use that to determine the initial waves, since we know each of the
% initial waves lies on the circular curve as the throat opens up to the
% max angle.

% Each wave is on the circular arc, spaced by THIS angle as the throat
% opens.
ddelta = delta_max/N;

% Instantiate arrays for the nodes
Kp_arr = zeros(N,N); % Left-running characteristics at nodes (invariants)
Km_arr = zeros(N,N); % Right-running characteristics at nodes (invariants)
M_arr = zeros(N,N); % Mach at nodes
delta_arr = zeros(N,N); % Flow Angles at nodes (deg)
nu_arr = zeros(N,N); % Prandtl-Meyer 'nu' at nodes (deg)
mu_arr = zeros(N,N); % Mach angles at nodes (deg)
x_arr = zeros(N,N); % x values of nodes (m)
y_arr = zeros(N,N); % y values of nodes (m)

%% Circular-Arc Throat Calculations (Task 1)
% Starting with the non-sharp throat comprised of a circular arc of radius
% equal to the throat radius. Discretize using the number of waves.

% Circular throat curve from horizontal to delta_max (degrees)
deltaThroatCircle(:, 1) = (0:ddelta:delta_max);

% Instantiate the 'nu' array with the same angle values
% True by definition of 'nu' since we start at M=1, so nu = delta on the
% throat expansion fan.
nuThroatCircle = deltaThroatCircle;

% And calcualte the mach angles upon which the characteristics go:
% first it needs to make mu a function of nu:
muFromNu = @(nu) asind(1./meyerMach(nu,gamma));

% Then calculate it.
% Remember it's useful because this determines angles of characteristic lines.
% mu = arcsin(1/M), and  dy/dx = tan(delta +/- mu)
muThroatCircle(:,1) = muFromNu(nuThroatCircle);

% Right-running (expansion waves that go down) characteristics from the
% throat (degrees)
% Definition of R: nu - delta = 0 for this initial expansion fan
KpThroatCircle = deltaThroatCircle + nuThroatCircle;
KmThroatCircle = nuThroatCircle - deltaThroatCircle; % should be zero

% Define an anonymous function to calculate the points
% of the circular-throat-opening from δ=0 to δ=δ_max

% Outputs: x, y in meters. Can take array or individual delta values.
% x,y coords along the circular throat opening to δ_max.
arcpoints = @(delta) deal( h_t .* sind(delta), h_t.*(2-cosd(delta)));
% CHECKED AGAINST DESMOS: https://www.desmos.com/calculator/x4hd9teecj

% Actually computes the coords of the circular part of the nozzle expanding.
[xThroatCircle, yThroatCircle] = arcpoints(deltaThroatCircle);

% Source points of the initial characteristic lines (waves)
nu0 = deltaThroatCircle(2:end);
theta0 = zeros(size(nu0));  % FIXED HIGH-MACH DIVERGENCE: start from sonic flow direction

% Originally was this but needed to delete due to high-mach divergence.
%theta0 = deltaThroatCircle(2:end);
%nu0 = theta0;

Kp0 = nu0 + theta0;
Km0 = nu0 - theta0;

x0 = xThroatCircle(2:end);
y0 = yThroatCircle(2:end);

mu0 = muThroatCircle(2:end);

%% Calculating Centerline Reflections (Task 1)
% The first big step to getting the nozzle geometry.

% For each initial throat wave 'i'
for i = 1:N
    % ----- THEORY -----
    % At centerline:
    % delta = 0
    % Kp is constant from source (derived in prior section)
    % So:
    %   nu = Kp
    %   Km = nu
    % Compute: M, mu
    % Compute slope
    % Intersect with y=0 (x-coordinate of first wave reflections)

    % Flow angle at centerline must zero due to symmetry
    delta_arr(i,1) = 0;

    % K+ was already calculated and stays the same
    Kp_arr(i,1) = Kp0(i); % From earlier array.
    
    % Then do the math: nu = (K+  + K- ) / 2
    % Where at centerline delta = 0
    % So K+  =  K-  = nu (reflection)
    nu_arr(i,1) = Kp_arr(i,1); % nu = K+
    Km_arr(i,1) = nu_arr(i,1); % K- = K+ = nu

    % Find Mach of the centerline given by 'nu'
    M_arr(i,1) = meyerMach(nu_arr(i,1), gamma);

    % Use that mach to find the mach angle
    mu_arr(i,1) = muFromNu(nu_arr(i,1));

    % Now we actually find the intersection since we know delta and mu
    % so we can find the slope using those:
    slope = tand(theta0(i) + mu0(i));
    dx = x0(i) + (-y0(i) / slope);
    
    if dx < x0(i)
        % flip to physically correct branch
        slope = tand(theta0(i) - mu0(i));
        dx = x0(i) + (-y0(i) / slope);
    end

    % And now find the intersection of the wave and centerline
    % Using the slope we just found
    % and save those coordinates:
    x_arr(i,1) = x0(i) + (-y0(i) / slope); % slope & initial y gives x to hit y=0.
    y_arr(i,1) = 0; % on centerline y=0.
end

%% Plot what we have so far (verify it looks ok)
debugGraphs = true;

if (debugGraphs)
    figure;
    hold on;
    grid on;
    axis equal;

    % Plot the circulr throat arc
    plot(xThroatCircle, yThroatCircle, 'b-', 'Linewidth', 2);

    % Plot the initial characteristic lines going from throat to center
    plot(x0, y0, 'ko', 'MarkerFaceColor', 'k');

    % Plot the centerline intersections
    plot(x_arr(:,1), y_arr(:,1), 'ro', 'MarkerFaceColor', 'r');

    % Plot the centerline itself:
    yline(0, 'k--');

    xlabel("x (m)");
    ylabel('y (m)');
    title('Check: Throat expansion to max angle, and first waves to centerline');

    legend('Throat arc', 'Initial wave sources', 'Centerline intersections', 'centerline');
end