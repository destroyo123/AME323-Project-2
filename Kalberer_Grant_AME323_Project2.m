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
% Requires having the a few toolboxes installed for certain functions

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



%% Setting up the arrays and variables
% todo: list what we put in here