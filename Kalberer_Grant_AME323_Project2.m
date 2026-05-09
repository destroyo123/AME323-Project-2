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
% * Prandtl-Meyer angle of local flow 'nu' (degrees)
% * Flow ratio of specific heats 'gamma,' usually 1.4
%
% 
%
% *Outputs:*
%
% * Local mach number 'M'
%
%
    meyerMach = @(nu) flowprandtlmeyer(gamma, nu, 'nu');
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
L_arr = zeros(N,N); % Left-running characteristics at nodes (invariants)
R_arr = zeros(N,N); % Right-running characteristics at nodes (invariants)
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
nuThroatCircle = deltaThroatCircle;

% And calcualte the mach angles upon which the characteristics go:
% first it needs to make mu a function of nu:
muFromNu = @(nu) asind(1/meyerMach(nu));

muThroatCircle = muFromNu(nuThroatCircle);

% Right-running (expansion waves that go down) characteristics from the
% throat (degrees)
RThroatCircle = deltaThroatCircle + nuThroatCircle;

% Define an anonymous function to calculate the points
% of the circular-throat-opening from δ=0 to δ=δ_max

% Outputs: x, y in meters. Can take array or individual delta values.
% x,y coords along the circular throat opening to δ_max.
arcpoints = @(delta) deal( h_t .* sind(delta), h_t.*(2-cosd(delta)));
% CHECKED AGAINST DESMOS: https://www.desmos.com/calculator/x4hd9teecj

% Actually computes the coords of the circular part of the nozzle expanding.
[xThroatCircle, yThroatCircle] = arcpoints(deltaThroatCircle);

% At the wall points along this arc:
R_arr(:,1) = RThroatCircle(2:length(RThroatCircle),1);
delta_arr(:,1) = deltaThroatCircle(2:length(RThroatCircle),1);
nu_arr(:,1) = delta_arr(:,1);
L_arr(:,1) = delta_arr(:,1) - nu_arr(:,1);

%% Calculating Centerline Reflections (Task 1)
% The first big step to getting the nozzle geometry.

% Throat can't be calculated exactly at M=1 or it breaks, so we give it a
% slight increase by 'eps'
eps = 1e-4;
M_arr(1,1) = 1+eps;
nu_arr(1,1) = 0; % deg
mu_arr(1,1) = 90; % deg

y_arr(1,1) = 0;
x_arr(1,1) = xThroatCircle(2,1) + ( y_arr(1,1) - yThroatCircle(2,1)  ) / tand((deltaThroatCircle(2,1) - muThroatCircle(2,1) - muThroatCircle(2,1))/2) ; 
x(1,1) = xarc(2,1) + (y(1,1) - yarc(2,1))/tand((ThetaArc(2,1) - MuArc(2,1) - MuArc(2,1))/2);

% Calculates the nodes on the first left (up) characteristic line "C+"
