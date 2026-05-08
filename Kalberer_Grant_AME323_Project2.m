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

% Given Relations

% Area ratio equation
AeAt = @(M) (1/M) * (2/(gamma+1) * (1 + (((gamma-1)/2) * M^2 )))^( (gamma+1)/(2*(gamma-1)) );

% Parameters

M_e = 6; % Exit mach
N_waves = 40; %Number of waves

h_n = 350; %Nozzle height in mm
h_n = h_n/1000; % Converted to m

h_t = (h_n ^2 / AeAt(M_e))^0.5; % Height of throat in meters



%% Setting up the arrays and contour stuff
% todo: list what we put in here

% Table for nozzle contour x, y, nu, delta, L, R or something


%% Functions
% Equations and Relations I mostly made for other homeworks & discussion
% problems. Equations come from NACA1135 unless otherwise stated.

%%
% *Prandtl-Meyer angle 'nu'*
% 
% *Inputs:*
%
%
% * Unperturbed Mach (AKA M1)
% * Ratio of specific heats, 'gamma' of flow, usually 1.4
% 
% *Outputs:*
%
% * Prandtl-Meyer angle 'nu' (degrees)
%
function prandtlMeyerAngle = nu(M, gamma)
    prandtlMeyerAngle = sqrt((gamma+1)/(gamma-1)) * atand(sqrt(((gamma-1)/(gamma+1))*(M.^2 - 1))) - atand(sqrt(M.^2 - 1));
    % OUTPUTS DEGREES
end
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
%
function mach = meyerMach(gamma, nu)
    mach = flowprandtlmeyer(gamma, nu, 'nu');
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