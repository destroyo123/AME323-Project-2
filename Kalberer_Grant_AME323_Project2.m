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

gray = [0, 0, 0];
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

N = 10; %Number of waves

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

% %% Calculating Centerline Reflections (Task 1)
% % The first big step to getting the nozzle geometry.
% 
% % For each initial throat wave 'i'
% for i = 1:N
%     % ----- THEORY -----
%     % At centerline:
%     % delta = 0
%     % Kp is constant from source (derived in prior section)
%     % So:
%     %   nu = Kp
%     %   Km = nu
%     % Compute: M, mu
%     % Compute slope
%     % Intersect with y=0 (x-coordinate of first wave reflections)
% 
%     % Flow angle at centerline must zero due to symmetry
%     delta_arr(i,1) = 0;
% 
%     % K+ was already calculated and stays the same
%     Kp_arr(i,1) = Kp0(i); % From earlier array.
% 
%     % Then do the math: nu = (K+  + K- ) / 2
%     % Where at centerline delta = 0
%     % So K+  =  K-  = nu (reflection)
%     nu_arr(i,1) = Kp_arr(i,1); % nu = K+
%     Km_arr(i,1) = nu_arr(i,1); % K- = K+ = nu
% 
%     % Find Mach of the centerline given by 'nu'
%     M_arr(i,1) = meyerMach(nu_arr(i,1), gamma);
% 
%     % Use that mach to find the mach angle
%     mu_arr(i,1) = muFromNu(nu_arr(i,1));
% 
%     % Now we actually find the intersection since we know delta and mu
%     % so we can find the slope using those:
%     slope = tand(theta0(i) + mu0(i));
%     dx = x0(i) + (-y0(i) / slope);
% 
%     if dx < x0(i)
%         % flip to physically correct branch
%         slope = tand(theta0(i) - mu0(i));
%         dx = x0(i) + (-y0(i) / slope);
%     end
% 
%     % And now find the intersection of the wave and centerline
%     % Using the slope we just found
%     % and save those coordinates:
%     x_arr(i,1) = x0(i) + (-y0(i) / slope); % slope & initial y gives x to hit y=0.
%     y_arr(i,1) = 0; % on centerline y=0.
% end

%% Calculating the initial waves
% This constructs the FIRST true MOC interaction row.
% We move from boundary "source points" to actual characteristic intersections.

% In MOC terms:
% (i,1) and (i+1,1) are boundary nodes from throat discretization
% We solve nodes (i,2) using intersections of:
%   C+ from (i,1)
%   C- from (i+1,1)

% BTW I caved and used AI for this because it would have teaken so damn
% long to rewrite myself and I still really don't get the nitty gritty..
% So this part onward is AI. (currently ChatGPT)

%% Initialize first interior column (j = 2)

for i = 1:N-1

    % ------------------------------------------------------------
    % STEP 1: Extract compatibility constants from neighbors
    % ------------------------------------------------------------

    % C+ comes from lower index (i)
    Kp = Kp0(i);

    % C- comes from upper index (i+1)
    Km = Km0(i+1);

    % ------------------------------------------------------------
    % STEP 2: Solve flow properties at intersection
    % ------------------------------------------------------------

    % Using your notation:
    % δ = (K+ - K-) / 2
    % ν = (K+ + K-) / 2

    delta_arr(i,2) = 0.5*(Kp - Km);
    nu_arr(i,2)    = 0.5*(Kp + Km);

    % Safety clamp (prevents PM function crash at high Mach)
    nu_arr(i,2) = max(min(nu_arr(i,2),130.45), 0);

    % Compute Mach and Mach angle
    M_arr(i,2)  = meyerMach(nu_arr(i,2), gamma);
    mu_arr(i,2) = muFromNu(nu_arr(i,2));

    % Store invariants at node
    Kp_arr(i,2) = Kp;
    Km_arr(i,2) = Km;

    % ------------------------------------------------------------
    % STEP 3: Compute geometry (intersection of characteristics)
    % ------------------------------------------------------------

    % C+ line from (i,1)
    m_plus = tand(delta_arr(i,1) + mu0(i));

    % C- line from (i+1,1)
    m_minus = tand(delta_arr(i+1,1) - mu0(i+1));

    % Points defining the two lines
    x1 = x0(i);     y1 = y0(i);
    x2 = x0(i+1);   y2 = y0(i+1);

    % Solve intersection
    x_int = (m_plus*x1 - m_minus*x2 + y2 - y1) / (m_plus - m_minus);
    y_int = m_plus*(x_int - x1) + y1;

    % Store node
    x_arr(i,2) = x_int;
    y_arr(i,2) = y_int;

end

%% CHECK YOURSELF!!
% Plot what we have so far (verify it looks ok)

debugGraphs = false;

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
    title('First Waves hitting centerline');

    legend('Throat arc', 'Initial wave sources', 'Centerline intersections', 'centerline');
end

%% ========================================================================
%  SECTION A: Full MOC Mesh Calculation
%  Centerline reflection, internal nodes, and wall points
%  ========================================================================

% -------------------------------------------------------------------------
% STEP 1: Propagate characteristics to the CENTERLINE (first reflection)
% -------------------------------------------------------------------------
% The C+ characteristics from the throat arc hit the centerline (y=0).
% At centerline: delta = 0, so nu = K+ (from the incoming C+ wave)
% After reflection, the wave becomes a C- wave with K- = nu (= K+_incoming)

% Arrays to store centerline nodes
x_center = zeros(N, 1);
y_center = zeros(N, 1);
nu_center = zeros(N, 1);
M_center = zeros(N, 1);
mu_center = zeros(N, 1);
Kp_center = zeros(N, 1);
Km_center = zeros(N, 1);
delta_center = zeros(N, 1);

for i = 1:N
    % At the centerline, delta = 0 (symmetry condition)
    delta_center(i) = 0;
    
    % K+ is preserved along C+ characteristic from source
    Kp_center(i) = Kp0(i);
    
    % At centerline: nu = K+ (since delta=0 => K+ = nu + delta = nu)
    nu_center(i) = Kp_center(i);
    
    % After reflection: K- = nu - delta = nu (since delta=0)
    Km_center(i) = nu_center(i);
    
    % Compute Mach number and Mach angle
    M_center(i) = meyerMach(nu_center(i), gamma);
    mu_center(i) = asind(1.0 / M_center(i));
    
    % Geometry: trace C+ line from source point (x0(i), y0(i)) to y = 0
    % Slope of C+ characteristic: dy/dx = tan(theta - mu) for a left-running
    % wave going DOWN toward centerline
    % The C+ characteristic slope heading toward centerline:
    slope_down = tand(theta0(i) - mu0(i));
    
    % If slope is positive or zero, use the other branch
    if slope_down >= 0
        slope_down = tand(theta0(i) + mu0(i));
        % Flip sign for going down
        slope_down = -abs(slope_down);
    end
    
    % Intersection with y = 0:
    % y0(i) + slope_down * (x_c - x0(i)) = 0
    x_center(i) = x0(i) + (-y0(i) / slope_down);
    y_center(i) = 0;
end

% -------------------------------------------------------------------------
% STEP 2: Build the full mesh using a structured node approach
% -------------------------------------------------------------------------
% We'll store ALL nodes in cell arrays for flexibility.
% The mesh is built column by column (each "column" = one set of 
% characteristic intersections moving from centerline toward wall).
%
% Structure:
%   - Start from centerline reflections (reflected C- waves go UP)
%   - These interact with the NEXT C+ wave coming from throat
%   - Continue until waves reach the wall
%   - At the wall: delta is chosen so that the reflected wave cancels
%     (wall condition: flow must be tangent to wall)

% We'll use a more general storage: node(row, col) where
%   col = 1: centerline points
%   col = 2..N: interior points moving toward wall
%   Last valid col for each row: wall point

% Maximum number of columns (reflections) we need
% For a minimum length nozzle: waves reflect off centerline then cancel at wall
% Total columns in mesh = N (centerline) + interior + wall
% For N initial waves, the mesh is triangular

% Let's build the full characteristic net
% We use a triangular mesh: 
%   Row i has (N - i + 1) nodes before hitting the wall
% Actually for MLN: each C+ from throat reflects off centerline, then 
% all reflected C- waves interact with each other and with subsequent C+ waves,
% finally terminating at the wall.

% I'll use a simpler approach: march through the net region by region.

% ---- Storage for ALL mesh nodes ----
maxNodes = N * (N + 1) / 2 + N + N; % upper bound
all_x = [];
all_y = [];
all_M = [];
all_nu = [];
all_delta = [];
all_mu = [];

% Also store characteristic LINE segments for plotting
% Each segment: [x1, y1, x2, y2, Mach_avg]
char_lines_Cplus = [];  % C+ (left-running) lines
char_lines_Cminus = []; % C- (right-running) lines

% ---- Rebuild mesh in a structured triangular grid ----
% node(i,j): i = row (which C+ family), j = column
% For MLN with N initial waves:
%   - There are N C+ characteristics from throat
%   - After centerline reflection, there are N C- characteristics
%   - Interior nodes form where C+ and C- cross

% Let's define:
%   Region 1 (Kernel): C+ from throat to centerline
%   Region 2 (Reflection/Transition): reflected C- interact with C+ 
%       and eventually define the wall

% =========================================================================
% KERNEL REGION: Throat to Centerline
% =========================================================================
% We already have:
%   Source points: (x0, y0) with flow properties (nu0, theta0, mu0, Kp0, Km0)
%   Centerline points: (x_center, y_center) with properties computed above

% Store source points as nodes
for i = 1:N
    all_x(end+1) = x0(i);
    all_y(end+1) = y0(i);
    all_M(end+1) = meyerMach(nu0(i), gamma);
    all_nu(end+1) = nu0(i);
    all_delta(end+1) = theta0(i);
    all_mu(end+1) = mu0(i);
end

% Store centerline points
for i = 1:N
    all_x(end+1) = x_center(i);
    all_y(end+1) = y_center(i);
    all_M(end+1) = M_center(i);
    all_nu(end+1) = nu_center(i);
    all_delta(end+1) = delta_center(i);
    all_mu(end+1) = mu_center(i);
end

% Store C+ lines from source to centerline
for i = 1:N
    M_avg = 0.5*(meyerMach(nu0(i), gamma) + M_center(i));
    char_lines_Cplus = [char_lines_Cplus; x0(i), y0(i), x_center(i), y_center(i), M_avg];
end

% =========================================================================
% INTERIOR NODES: Between initial C+ lines (already partially computed)
% =========================================================================
% We already computed column j=2 nodes in your code above.
% Let's also store those and add C+ line segments from them.

% Store the j=2 interior nodes (between adjacent C+ waves, before centerline)
for i = 1:N-1
    if x_arr(i,2) > 0 && y_arr(i,2) > 0
        all_x(end+1) = x_arr(i,2);
        all_y(end+1) = y_arr(i,2);
        all_M(end+1) = M_arr(i,2);
        all_nu(end+1) = nu_arr(i,2);
        all_delta(end+1) = delta_arr(i,2);
        all_mu(end+1) = mu_arr(i,2);
    end
end

% =========================================================================
% REFLECTION REGION: After centerline, C- waves go up toward wall
% =========================================================================
% After reflecting off centerline, each wave i becomes a C- characteristic
% with K- = Km_center(i) going upward.
% These C- waves interact with other C+ waves (from later throat emissions)
% and with each other.

% For a MINIMUM LENGTH NOZZLE, the wall is shaped so that ALL reflected 
% waves are cancelled at the wall (no secondary reflection off wall needed).
% This means each C- from centerline goes directly to the wall and determines
% the wall angle there.

% However, you asked for at least one centerline reflection. Let's do:
%   1) C+ from throat -> hits centerline (done above)
%   2) Reflected C- from centerline -> interacts with other C+ in interior
%   3) C- waves reach wall -> wall angle cancels them (MLN condition)

% ---- Compute interior interaction nodes in the reflection region ----
% After centerline reflection, we have N upgoing C- waves.
% These interact with C+ waves that are still propagating.

% For the triangular net after reflection:
% Node(i,j) in reflection region:
%   C- comes from centerline reflection point i (going up-right)
%   C+ comes from centerline reflection point j (going up-right... no)
% 
% Actually, in a standard MLN MOC net:
%   After the kernel (throat to centerline), the "transition" region has
%   nodes where reflected C- waves cross the original C+ waves.

% Let me restructure using a cleaner indexing:
% 
% The net after centerline has these nodes:
%   For reflected C- wave from center point i (i=1..N):
%     It crosses C+ waves from throat points i+1, i+2, ..., N
%     Then it hits the wall.
%
% Node(i,j) in transition region:
%   i = index of the C- wave (from center reflection i)  
%   j = index of the C+ wave it's crossing

% Transition region nodes
% x_trans(i,j), y_trans(i,j): node where C-(from center i) crosses C+(from throat j)
% Valid for j > i (C- from earlier reflection crosses later C+ waves)

x_trans = NaN(N, N);
y_trans = NaN(N, N);
nu_trans = NaN(N, N);
delta_trans = NaN(N, N);
M_trans = NaN(N, N);
mu_trans = NaN(N, N);
Kp_trans = NaN(N, N);
Km_trans = NaN(N, N);

% We need to march through the net carefully.
% The approach: process nodes column by column (by C+ index j)
% For each j, process all C- waves (i < j) that cross it.

% But first, let's establish what we know at the boundaries:
% 
% Along each C+ (index j, from throat):
%   K+ = Kp0(j) = nu0(j) + theta0(j)  [constant along this C+]
%
% Along each C- (index i, from centerline reflection):
%   K- = Km_center(i) = nu_center(i)   [constant along this C-]

% Interior node where C+(j) meets C-(i), with i < j:
%   nu = (K+ + K-) / 2 = (Kp0(j) + Km_center(i)) / 2
%   delta = (K+ - K-) / 2 = (Kp0(j) - Km_center(i)) / 2

% For geometry, we need to know the PREVIOUS node on each characteristic:
%   - Previous node on C+(j): 
%       if i == j-1: it's the centerline point j (x_center(j), y_center(j))
%       else: it's the transition node (i-1, j)
%   - Previous node on C-(i):
%       if j == i+1: it's the centerline point i (x_center(i), y_center(i))  
%       else: it's the transition node (i, j-1)

% March through the net
for j = 2:N       % C+ wave index (from throat)
    for i = 1:j-1  % C- wave index (from centerline reflection)
        
        % Flow properties from compatibility equations
        Kp_trans(i,j) = Kp0(j);           % K+ constant along C+(j)
        Km_trans(i,j) = Km_center(i);     % K- constant along C-(i)
        
        nu_trans(i,j) = 0.5 * (Kp_trans(i,j) + Km_trans(i,j));
        delta_trans(i,j) = 0.5 * (Kp_trans(i,j) - Km_trans(i,j));
        
        % Clamp nu for safety
        nu_trans(i,j) = max(min(nu_trans(i,j), 130.45), 0.01);
        
        M_trans(i,j) = meyerMach(nu_trans(i,j), gamma);
        mu_trans(i,j) = asind(1.0 / M_trans(i,j));
        
        % ---- Geometry: find (x,y) of this node ----
        % Previous node on C+(j):
        if i == 1
            % First crossing of this C+ after centerline
            % Previous point is centerline node j
            xA = x_center(j);
            yA = y_center(j);
            delta_A = delta_center(j);
            mu_A = mu_center(j);
        else
            % Previous point is transition node (i-1, j)
            xA = x_trans(i-1, j);
            yA = y_trans(i-1, j);
            delta_A = delta_trans(i-1, j);
            mu_A = mu_trans(i-1, j);
        end
        
        % Previous node on C-(i):
        if j == i + 1
            % First crossing of this C- after centerline
            % Previous point is centerline node i
            xB = x_center(i);
            yB = y_center(i);
            delta_B = delta_center(i);
            mu_B = mu_center(i);
        else
            % Previous point is transition node (i, j-1)
            xB = x_trans(i, j-1);
            yB = y_trans(i, j-1);
            delta_B = delta_trans(i, j-1);
            mu_B = mu_trans(i, j-1);
        end
        
        % C+ slope (from point A going up): dy/dx = tan(delta + mu)
        % Actually for C+ (left-running), slope = tan(delta - mu) 
        % (it runs from lower-left to upper-right in supersonic flow)
        % Convention: C+ has slope tan(theta - mu), C- has slope tan(theta + mu)
        % But since C+ is going from centerline UP, and C- from center UP too...
        
        % Let me be careful with signs:
        % C+ characteristic (left-running): slope = tan(delta - mu)  [goes up-right from centerline]
        % C- characteristic (right-running): slope = tan(delta + mu) [goes up-right from centerline]
        
        % Average the slopes (predictor-corrector style):
        % For C+ from A:
        slope_Cplus_A = tand(delta_A - mu_A);
        slope_Cplus_node = tand(delta_trans(i,j) - mu_trans(i,j));
        slope_Cplus = 0.5*(slope_Cplus_A + slope_Cplus_node);
        
        % For C- from B:
        slope_Cminus_B = tand(delta_B + mu_B);
        slope_Cminus_node = tand(delta_trans(i,j) + mu_trans(i,j));
        slope_Cminus = 0.5*(slope_Cminus_B + slope_Cminus_node);
        
        % Solve intersection:
        % yA + slope_Cplus*(x - xA) = yB + slope_Cminus*(x - xB)
        if abs(slope_Cplus - slope_Cminus) < 1e-12
            % Parallel lines - shouldn't happen but safeguard
            x_trans(i,j) = 0.5*(xA + xB);
            y_trans(i,j) = 0.5*(yA + yB);
        else
            x_trans(i,j) = (yB - yA + slope_Cplus*xA - slope_Cminus*xB) / (slope_Cplus - slope_Cminus);
            y_trans(i,j) = yA + slope_Cplus*(x_trans(i,j) - xA);
        end
        
        % Store for plotting
        all_x(end+1) = x_trans(i,j);
        all_y(end+1) = y_trans(i,j);
        all_M(end+1) = M_trans(i,j);
        all_nu(end+1) = nu_trans(i,j);
        all_delta(end+1) = delta_trans(i,j);
        all_mu(end+1) = mu_trans(i,j);
        
        % Store characteristic line segments
        M_avg_cp = 0.5*(meyerMach(max(0.01, 0.5*(Kp_trans(i,j)+Km_trans(i,j))), gamma) + ...
                        M_trans(i,j));
        char_lines_Cplus = [char_lines_Cplus; xA, yA, x_trans(i,j), y_trans(i,j), M_avg_cp];
        char_lines_Cminus = [char_lines_Cminus; xB, yB, x_trans(i,j), y_trans(i,j), M_avg_cp];
    end
end


% =========================================================================
% WALL POINTS: Where each C- wave terminates (cancellation condition)
% =========================================================================
% At the wall, the flow must be tangent to the wall.
% For a MINIMUM LENGTH NOZZLE, the wall angle equals the local flow angle.
% The wall CANCELS the wave: no reflection back.
%
% For C-(i) hitting the wall:
%   The last interior node before the wall on C-(i) is:
%     If i < N: node (i, N) in transition region
%     If i = N: this C- only comes from centerline and goes straight to wall
%               (no interior crossings since there's no C+ with index > N)
%
% At wall point for C-(i):
%   K- = Km_center(i) is still constant
%   Wall condition: delta_wall = delta (flow tangent to wall)
%   For cancellation (no reflected wave): 
%     The wall angle equals the flow deflection angle at that point.
%   
%   From the LAST C+ that this C- crossed (which is C+(N) for i < N,
%   or no crossing for i = N):
%     K+ = Kp0(N) for nodes where j=N was the last crossing
%   
%   For i = N (the last reflected wave): it goes directly from centerline to wall
%     K+ must be determined by the wall condition.
%     At wall: delta_wall = delta = (K+ - K-)/2
%     For MLN: the wall cancels the wave, meaning:
%       The wall slope = delta at that point
%       K+ at wall = K+(from the last C+ it would have crossed)
%       But for the Nth wave, it hasn't crossed any C+ after reflection.
%       So we use: at wall, delta = nu (cancellation means Km = 0 effectively)
%       Actually: For MLN, wall angle = delta = (Kp - Km)/2
%       and nu = (Kp + Km)/2
%       With wall cancellation: the reflected K+ would equal K- 
%       => no new wave. This means wall angle = delta at that point.

% Wall point arrays
x_wall = zeros(N, 1);
y_wall = zeros(N, 1);
delta_wall = zeros(N, 1);
M_wall = zeros(N, 1);
nu_wall = zeros(N, 1);
mu_wall = zeros(N, 1);

for i = 1:N
    % The C- wave from centerline reflection i
    % K- is constant = Km_center(i)
    Km_w = Km_center(i);
    
    % For MLN wall cancellation:
    % At wall, the condition is that the wall absorbs the wave.
    % This means: K+ (of the would-be reflected wave) = K- 
    % => delta_wall = 0?? No...
    %
    % Actually for MLN: 
    % The exit flow is uniform at M_e with delta=0.
    % The LAST characteristic to reach the wall defines the exit.
    % The wall angle at each point = flow angle delta at that point.
    % 
    % For wall point of C-(i):
    %   It has crossed all C+ waves up to C+(N).
    %   After crossing C+(N), it goes to the wall.
    %   At that last interior node (i,N): we know nu, delta, etc.
    %   The wall point: K+ = Kp0_wall (from wall compatibility)
    %   Wall condition: delta_wall = flow angle at wall
    %   For cancellation: nu_wall = K+ at wall (like centerline but for wall)
    %   => K+_wall = nu_wall + delta_wall
    %   And K- = nu_wall - delta_wall = Km_center(i)
    %   For the wave to be cancelled: no new K+ is generated
    %   This means: wall slope = delta at the wall point
    %   And nu_wall and delta_wall come from the last interior node's C+ value.
    
    % Simpler approach for MLN:
    % All waves reach the wall after the transition region.
    % At the wall: K+ = nu_e + 0 = nu_e (exit condition: uniform flow)
    % Wait no. Let me think again...
    
    % For a minimum length nozzle, the EXIT has uniform flow at M_e, delta=0.
    % The wall is shaped such that each C- wave, upon hitting the wall,
    % is perfectly cancelled (turned to zero deflection by the exit).
    % 
    % The K+ value along the wall changes from point to point.
    % At each wall point: 
    %   K- = Km_center(i) (constant along C-)
    %   K+ comes from the last C+ characteristic that the C- crossed.
    %   For i <= N-1: the last C+ crossed is C+(N), so K+ = Kp0(N)
    %   For i = N: it hasn't crossed any C+ after reflection, 
    %              the wall condition gives K+ = nu_e (exit PM angle)
    %              Actually K+ = Kp0(N) = nu0(N) + theta0(N) for all
    
    % Let me use: at wall point for C-(i):
    %   If i < N: last interior node is (i, N) in transition region
    %   If i = N: last node is centerline point N
    
    % Actually the simplest correct approach:
    % The wall point for C-(i) after crossing all available C+ waves:
    % K+ at wall = Kp0(N) for all i (since all C- waves eventually 
    % cross C+(N), the last and strongest C+ wave)
    % EXCEPT for i=N which reflects from centerline and goes directly to wall
    % without crossing any C+ (since there's no C+(N+1))
    
    % For i = N: 
    %   K- = Km_center(N) = nu_center(N) = Kp0(N) = nu0(N) + theta0(N)
    %   At wall: we need K+ from somewhere. Since this is the LAST wave,
    %   at the exit the flow should be uniform: delta=0, nu=nu_e
    %   So K+ = nu_e + 0 = nu_e, K- = nu_e - 0 = nu_e
    %   Indeed Km_center(N) should equal nu_e for MLN!
    %   (Because Kp0(N) = nu0(N) + theta0(N) = delta_max + 0 = delta_max
    %    Wait: theta0 = 0 (you set it to zero), nu0 = delta_max for last wave
    %    So Kp0(N) = delta_max + 0 = delta_max
    %    And nu_center(N) = Kp0(N) = delta_max = nu_e/2
    %    Hmm, that's nu_e/2, not nu_e...)
    
    % OK let me just use the straightforward MOC approach:
    % For each wall point, the flow properties are determined by
    % the incoming C- characteristic, and the wall condition that 
    % the flow is tangent to the wall (delta = wall angle).
    % For MLN: at the wall, the C- is cancelled, meaning:
    %   The wall turns the flow so that no reflected C+ is generated.
    %   Condition: delta_wall is such that K+_reflected = K+_incoming
    %   In practice for MLN: delta_wall = delta at that node.
    
    % Let's just use the last known node on each C- and trace to wall.
    % The "wall" K+ equals the K+ of the last C+ wave that was crossed.
    
    if i < N
        % Last interior node on C-(i) is (i, N)
        % After that, it goes to the wall
        % K+ at wall = Kp0(N) (still on the last C+ family... 
        % Actually no: after the transition node (i,N), the C- continues
        % to the wall. The wall point has:
        %   K- = Km_center(i) [still constant]
        %   K+ = Kp0(N) [from the last C+ crossed at node (i,N)]
        % So flow at wall:
        Kp_w = Kp0(N);
        nu_wall(i) = 0.5*(Kp_w + Km_w);
        delta_wall(i) = 0.5*(Kp_w - Km_w);
        
        % Previous node on C-(i) going to wall:
        xB = x_trans(i, N);
        yB = y_trans(i, N);
        delta_B = delta_trans(i, N);
        mu_B = mu_trans(i, N);
    else
        % i = N: goes directly from centerline to wall (no C+ crossings after)
        % For the last wave in MLN, at the wall the flow becomes the exit flow
        % K- = Km_center(N), and the wall condition with uniform exit:
        % Actually, since this is the last wave and defines the exit,
        % at this wall point: nu = nu_e, delta = 0 (exit condition)
        % Check: K- = nu - delta = nu_e, K+ = nu + delta = nu_e
        % But Km_center(N) = nu_center(N) = Kp0(N) = nu0(N) = delta_max = nu_e/2
        % This doesn't match nu_e unless we have more reflections...
        
        % For a SINGLE reflection MLN (which is what we're doing):
        % The exit Mach is determined by the characteristics, not prescribed.
        % With theta0=0 and nu0 = delta values from 1*ddelta to N*ddelta:
        % Kp0(N) = nu0(N) + theta0(N) = N*ddelta = delta_max = nu_e/2
        % nu_center(N) = Kp0(N) = nu_e/2
        % So the ACTUAL exit Mach from this single-reflection net is determined
        % by nu = nu_e/2... unless we account for multiple reflections.
        
        % For a TRUE MLN: delta_max = nu_e/2 and the exit nu = nu_e.
        % The wall must turn the flow back to delta=0 while nu increases to nu_e.
        % This happens because at the wall point:
        %   K+ (from the wall) = nu + delta, and since the wall cancels,
        %   the K+ leaving the wall = K+ arriving (no new wave).
        %   But the wall angle IS delta, so the flow follows the wall.
        
        % For the last C- (i=N):
        %   At centerline: nu = nu_e/2, delta = 0
        %   At wall (cancellation): delta_wall = 0 (exit), nu_wall = nu_e/2
        %   Hmm, that gives exit Mach from nu_e/2, not nu_e...
        
        % I think the issue is that with a single centerline reflection and
        % N initial waves, the exit Mach = meyerMach(nu_e/2 + delta_max/something)
        % Let me reconsider: 
        % For the wall points in standard MLN:
        %   At wall point where C-(i) arrives:
        %     K- = Km_center(i) = nu_center(i)
        %     Wall cancellation means: delta_wall = (K+_last - K-)/2
        %     And the wall point gets its K+ from the PREVIOUS wall point's
        %     flow condition... This is actually an iterative wall construction.
        
        % SIMPLER CORRECT APPROACH:
        % In a standard 2D planar MLN with the sharp-corner throat approx:
        %   - Initial expansion fan from theta=0 to theta=theta_max at throat
        %   - theta_max = nu_e/2 (this is why it's called "half the PM angle")
        %   - Waves hit centerline and reflect
        %   - Reflected waves hit the wall where they are cancelled
        %   - At each wall point: 
        %       K- is known (from reflection)
        %       Wall condition: flow tangent to wall => delta = wall slope
        %       Cancellation: K+ = K+ from the neighboring wall point
        %       For the FIRST wall point (i=N, farthest downstream on C+ line N):
        %         K+ comes from the initial expansion: Kp0(N)
        
        % For i=N: previous node is centerline point N
        Kp_w = Kp0(N); % Since the last C+ is from the throat
        % Actually for i=N from centerline, there's no C+ to cross,
        % but the K+ at this node was Kp0(N) = Kp_center(N) anyway.
        nu_wall(i) = 0.5*(Kp_w + Km_w);
        delta_wall(i) = 0.5*(Kp_w - Km_w);
        
        xB = x_center(N);
        yB = y_center(N);
        delta_B = delta_center(N);
        mu_B = mu_center(N);
    end
    
    % Clamp
    nu_wall(i) = max(min(nu_wall(i), 130.45), 0.01);
    M_wall(i) = meyerMach(nu_wall(i), gamma);
    mu_wall(i) = asind(1.0 / M_wall(i));
    
    % Trace C- from previous node to wall point
    % C- slope: tan(delta + mu)
    slope_Cminus_B = tand(delta_B + mu_B);
    slope_Cminus_wall = tand(delta_wall(i) + mu_wall(i));
    slope_Cm = 0.5*(slope_Cminus_B + slope_Cminus_wall);
    
    % Wall point lies along this C- line from (xB, yB).
    % We also need the wall curve. For the FIRST wall point, 
    % we can determine it from the previous wall point or from geometry.
    % For MLN: the wall starts at the end of the throat arc 
    % (at delta_max) and each subsequent wall point is determined
    % by the C- intersection with a line of slope tan(delta_wall(i)) from
    % the previous wall point.
    
    % For the first wall point (i=N is actually the closest to throat):
    % Let me reorder: wall point from C-(1) is farthest downstream,
    % wall point from C-(N) is closest to throat.
    % Actually C-(1) reflects first (smallest nu), C-(N) reflects last (largest nu).
    % C-(N) is the strongest wave and hits the wall closest to throat.
    % C-(1) is the weakest and hits the wall farthest from throat.
    
    % For now, let's just trace each C- to find the wall y-coordinate
    % by matching with the wall growth condition.
    % SIMPLE APPROACH: assume the wall point is where the C- line from 
    % the last interior node intersects with a line from the previous 
    % wall point at the wall angle.
    
    % Store wall point (we'll compute geometry after all properties are known)
    % For now, trace the C- and assume wall is "far enough" 
    % We'll intersect with wall lines next.
end

% =========================================================================
% WALL GEOMETRY CONSTRUCTION
% =========================================================================
% Build wall from throat arc endpoint to exit.
% Process wall points from closest to throat (i=N) to farthest (i=1).
% Wall starts at the end of the circular arc: (xThroatCircle(end), yThroatCircle(end))
% with slope = tan(delta_max).

% Reorder: wall point N is closest to throat, wall point 1 is farthest
% (C-(N) has the shortest path after reflection)
% Actually let me re-examine: C-(i) from center point i.
% Center point 1 has the smallest nu (weakest), center point N has largest.
% After reflection, C-(1) goes up with small nu (almost sonic),
% C-(N) goes up with large nu (high Mach).
% In terms of x-position: center point 1 is closest to throat,
% center point N is farthest downstream.
% So C-(1) reflects closest to throat and hits wall closest to throat,
% C-(N) reflects farthest and hits wall farthest downstream.

% Wall starting point (end of circular arc)
x_wall_start = xThroatCircle(end);
y_wall_start = yThroatCircle(end);
delta_wall_start = delta_max;

% Build wall points sequentially
x_wall_pts = zeros(N+1, 1);
y_wall_pts = zeros(N+1, 1);
delta_wall_pts = zeros(N+1, 1);

x_wall_pts(1) = x_wall_start;
y_wall_pts(1) = y_wall_start;
delta_wall_pts(1) = delta_wall_start;

for i = 1:N
    % C-(i) hits the wall at point i+1
    % Previous node on C-(i) before wall:
    if i < N && ~isnan(x_trans(i, N))
        xB = x_trans(i, N);
        yB = y_trans(i, N);
        delta_B = delta_trans(i, N);
        mu_B = mu_trans(i, N);
    elseif i < N
        % Fallback: use last valid transition node
        last_valid_j = find(~isnan(x_trans(i,:)), 1, 'last');
        if ~isempty(last_valid_j)
            xB = x_trans(i, last_valid_j);
            yB = y_trans(i, last_valid_j);
            delta_B = delta_trans(i, last_valid_j);
            mu_B = mu_trans(i, last_valid_j);
        else
            xB = x_center(i);
            yB = y_center(i);
            delta_B = delta_center(i);
            mu_B = mu_center(i);
        end
    else
        % i = N: comes directly from centerline
        xB = x_center(N);
        yB = y_center(N);
        delta_B = delta_center(N);
        mu_B = mu_center(N);
    end
    
    % Slope of C- arriving at wall: average of departure and arrival slopes
    slope_Cm_B = tand(delta_B + mu_B);
    slope_Cm_W = tand(delta_wall(i) + mu_wall(i));
    slope_Cm = 0.5*(slope_Cm_B + slope_Cm_W);
    
    % Wall line from previous wall point with slope = tan(delta_wall(i))
    % (The wall angle at the NEW point is delta_wall(i))
    % Use average wall angle between previous and current point:
    slope_wall = tand(0.5*(delta_wall_pts(i) + delta_wall(i)));
    
    % Intersection of:
    %   C- line: y = yB + slope_Cm * (x - xB)
    %   Wall line: y = y_wall_pts(i) + slope_wall * (x - x_wall_pts(i))
    
    if abs(slope_Cm - slope_wall) < 1e-12
        % Nearly parallel - advance slightly
        x_wall_pts(i+1) = x_wall_pts(i) + 0.01*h_t;
        y_wall_pts(i+1) = y_wall_pts(i) + slope_wall*0.01*h_t;
    else
        x_wall_pts(i+1) = (yB - y_wall_pts(i) + slope_wall*x_wall_pts(i) - slope_Cm*xB) / ...
                           (slope_wall - slope_Cm);
        y_wall_pts(i+1) = yB + slope_Cm*(x_wall_pts(i+1) - xB);
    end
    
    delta_wall_pts(i+1) = delta_wall(i);
    
    % Store characteristic line to wall
    M_avg = 0.5*(meyerMach(max(0.01,nu_wall(i)),gamma) + M_wall(i));
    char_lines_Cminus = [char_lines_Cminus; xB, yB, x_wall_pts(i+1), y_wall_pts(i+1), M_avg];
    
    % Store wall node
    all_x(end+1) = x_wall_pts(i+1);
    all_y(end+1) = y_wall_pts(i+1);
    all_M(end+1) = M_wall(i);
    all_nu(end+1) = nu_wall(i);
    all_delta(end+1) = delta_wall(i);
    all_mu(end+1) = mu_wall(i);
end

% Add final straight section to exit (uniform flow, delta=0)
% The last wall point should have delta ≈ 0. If not, extend to where delta=0.
if abs(delta_wall_pts(end)) > 0.1
    % Add exit point where flow becomes uniform
    % Extend wall straight at last delta until reaching exit height
    x_exit = x_wall_pts(end) + 2*h_t; % arbitrary extension
    y_exit = y_wall_pts(end) + tand(delta_wall_pts(end))*(x_exit - x_wall_pts(end));
    x_wall_pts(end+1) = x_exit;
    y_wall_pts(end+1) = y_exit;
    delta_wall_pts(end+1) = 0;
end

% =========================================================================
% BUILD CONNECTED WAVE POLYLINES (for clean plotting, no zigzag)
% =========================================================================
% Each C+ wave (index j): throat source j -> centerline j -> up through
%   transition nodes (i,j) for i=1..j-1 -> wall point j
% Each C- wave (index i): centerline i -> transition nodes (i,j) for j=i+1..N -> wall point i+1

% C+ polylines: one cell per wave j=1..N
Cplus_poly = cell(N, 1);
for j = 1:N
    xs = x0(j);   ys = y0(j);   % start: throat source
    xc = x_center(j); yc = y_center(j); % centerline
    % Collect transition nodes on this C+ (i goes from 1 to j-1, but
    % in transition region C+(j) is crossed by C-(i) for i < j.
    % The nodes on C+(j) going upward are (1,j), (2,j),...,(j-1,j))
    x_path = [xs; xc];
    y_path = [ys; yc];
    for i = 1:j-1
        if ~isnan(x_trans(i,j))
            x_path(end+1) = x_trans(i,j);
            y_path(end+1) = y_trans(i,j);
        end
    end
    % Wall point for C+(j): the C-(j) from centerline hits wall at x_wall_pts(j+1)
    % Actually C+(j) terminates at the wall where C-(j) hits it: wall point index j+1
    % Wait - C+(j) goes throat->centerline(j). After that it's "done" as a C+.
    % The upward segments above centerline are NOT C+(j), they're C-(i) crossing C+(j).
    % So C+(j) polyline is just: source(j) -> centerline(j)
    Cplus_poly{j} = [x0(j), y0(j); x_center(j), y_center(j)];
end

% C- polylines: one cell per wave i=1..N
% C-(i) goes: centerline(i) -> trans(i,i+1) -> trans(i,i+2) -> ... -> trans(i,N) -> wall(i+1)
Cminus_poly = cell(N, 1);
for i = 1:N
    x_path = x_center(i);
    y_path = y_center(i);
    for j = i+1:N
        if ~isnan(x_trans(i,j))
            x_path(end+1) = x_trans(i,j);
            y_path(end+1) = y_trans(i,j);
        end
    end
    % Terminate at wall point i+1
    if i+1 <= length(x_wall_pts)
        x_path(end+1) = x_wall_pts(i+1);
        y_path(end+1) = y_wall_pts(i+1);
    end
    Cminus_poly{i} = [x_path(:), y_path(:)];
end

fprintf('--- MOC Mesh Calculation Complete ---\n');
fprintf('Number of initial waves: %d\n', N);
fprintf('Exit Mach (design): %.2f\n', M_e);
fprintf('Max wall angle (delta_max): %.2f deg\n', delta_max);
fprintf('Throat half-height: %.4f m\n', h_t);
fprintf('Number of wall points: %d\n', length(x_wall_pts));
fprintf('Number of characteristic segments (C+): %d\n', size(char_lines_Cplus,1));
fprintf('Number of characteristic segments (C-): %d\n', size(char_lines_Cminus,1));

%% ========================================================================
%  SECTION B: Plot Nozzle with Patch-Based Mach Fill + Clean Characteristic Lines
%  ========================================================================

figure('Name', 'MOC Minimum Length Nozzle', 'Position', [100 100 1400 600]);
hold on;
grid on;
axis equal;

M_min_plot = 1.0;
M_max_plot = max(all_M(:));
if M_max_plot <= M_min_plot; M_max_plot = M_min_plot + 1; end

colormap('jet');
caxis([M_min_plot M_max_plot]);

% -------------------------------------------------------------------------
% PATCH-BASED MACH FILL
% -------------------------------------------------------------------------
% Strategy: fill the nozzle using patch() triangles/quads between
% characteristic nodes. Each patch gets the average Mach of its corners.
% This is exact per-region coloring with no interpolation artifacts.
%
% Regions to fill:
%  1) Kernel triangles: between adjacent C+ waves, from throat to centerline
%  2) Transition quads: between C- waves crossing C+ waves
%  3) Wall triangles: between last transition node, wall point, and neighbor

n_cmap = 256;
cmap = colormap;  % 64 colors by default; force 256
colormap(jet(n_cmap));
cmap = colormap;

getColor = @(M_val) cmap(max(1, min(n_cmap, ...
    round((M_val - M_min_plot)/(M_max_plot - M_min_plot)*(n_cmap-1)) + 1)), :);

% -------------------------------------------------------------------------
% PATCH FILL: color each quadrilateral cell between characteristic waves
% -------------------------------------------------------------------------
% The MOC mesh forms a triangular grid. We fill it by identifying all 
% quadrilateral (or triangular) cells and patch-filling each one.
%
% The nodes are:
%   Source points:     (x0(j),      y0(j))       for j=1..N  [on throat arc]
%   Centerline points: (x_center(i), y_center(i)) for i=1..N  [on y=0]
%   Transition nodes:  (x_trans(i,j), y_trans(i,j)) for j>i   [interior]
%   Wall points:       (x_wall_pts(k), y_wall_pts(k)) for k=1..N+1
%
% Cell types:
%   A) Kernel quads: bounded by C+(j), C+(j+1), and the centerline.
%      Corners: source(j), source(j+1), centerline(j+1), centerline(j)
%   B) Transition quads: bounded by adjacent C- waves crossing adjacent C+ waves.
%      Corners: trans(i,j), trans(i+1,j), trans(i+1,j+1), trans(i,j+1)
%      (with boundary substitutions for i=0: centerline, and j=N: wall)
%   C) Left-edge triangles: centerline(i) to trans(i,i+1) to centerline(i+1)
%   D) Wall-edge quads: trans(i,N) to wall(i+1) to wall(i+2) to trans(i+1,N)

% Helper: get Mach at a node type
getMach_src    = @(j) meyerMach(nu0(j), gamma);
getMach_center = @(i) M_center(i);
getMach_trans  = @(i,j) M_trans(i,j);
getMach_wall   = @(k) M_wall(min(k, N));  % wall point k+1 has M_wall(k)

fillQuad = @(x4, y4, M4) patch(x4, y4, mean(M4), ...
    'EdgeColor', 'none', 'FaceColor', getColor(mean(M4)));
fillTri  = @(x3, y3, M3) patch(x3, y3, mean(M3), ...
    'EdgeColor', 'none', 'FaceColor', getColor(mean(M3)));

% ---- A) Kernel region: quads between adjacent C+ waves ----
% Each quad: source(j), source(j+1), centerline(j+1), centerline(j)
% (Note: centerline points are ordered by x, so centerline(1) is closest
%  to throat and centerline(N) is farthest. The quad sweeps the full kernel.)
for j = 1:N-1
    x4 = [x0(j),        x0(j+1),        x_center(j+1),  x_center(j)];
    y4 = [y0(j),        y0(j+1),        y_center(j+1),  y_center(j)];
    M4 = [getMach_src(j), getMach_src(j+1), getMach_center(j+1), getMach_center(j)];
    fillQuad(x4, y4, M4);
end

% ---- B) Transition region: quads between C- waves ----
% Node indexing: row=i (C- index), col=j (C+ index crossed), j>i
% For a full quad at (i,j):
%   Bottom-left:  if i==1  -> centerline(j),   else trans(i-1, j)
%   Bottom-right: if i==1  -> centerline(j+1), else trans(i-1, j+1)  [if j<N]
%   Top-right:    trans(i, j+1)  [if j<N], else wall point
%   Top-left:     trans(i, j)

% We loop over all valid quads where both (i,j) and (i,j+1) exist
for i = 1:N-1
    for j = i+1:N-1
        if isnan(x_trans(i,j)) || isnan(x_trans(i,j+1)); continue; end
        
        % Top edge of quad: trans(i,j) and trans(i,j+1)
        xTL = x_trans(i,j);   yTL = y_trans(i,j);   MTL = M_trans(i,j);
        xTR = x_trans(i,j+1); yTR = y_trans(i,j+1); MTR = M_trans(i,j+1);
        
        % Bottom edge: trans(i-1,j) and trans(i-1,j+1), or centerline if i==1
        % Wait - in this mesh "up" means toward wall, so:
        % trans(i,j) is ABOVE trans(i-1,j) (higher i = higher up/more reflected)
        % Re-check: i is C- index (reflected from centerline(i)),
        %           j is C+ index crossed. Higher i = farther from throat on centerline
        %           = lower y position of node? Let's just use the 4 corners directly.
        
        if i == 1
            xBL = x_center(j);   yBL = y_center(j);   MBL = M_center(j);
            xBR = x_center(j+1); yBR = y_center(j+1); MBR = M_center(j+1);
        else
            if isnan(x_trans(i-1,j)) || isnan(x_trans(i-1,j+1)); continue; end
            xBL = x_trans(i-1,j);   yBL = y_trans(i-1,j);   MBL = M_trans(i-1,j);
            xBR = x_trans(i-1,j+1); yBR = y_trans(i-1,j+1); MBR = M_trans(i-1,j+1);
        end
        
        x4 = [xBL, xBR, xTR, xTL];
        y4 = [yBL, yBR, yTR, yTL];
        M4 = [MBL, MBR, MTR, MTL];
        fillQuad(x4, y4, M4);
    end
end

% ---- C) Left-edge triangles of transition region ----
% These are the triangles between centerline(i), centerline(i+1), and trans(i, i+1)
% (the very first cell each C- creates after reflecting)
for i = 1:N-1
    if isnan(x_trans(i, i+1)); continue; end
    x3 = [x_center(i), x_center(i+1), x_trans(i, i+1)];
    y3 = [y_center(i), y_center(i+1), y_trans(i, i+1)];
    M3 = [M_center(i), M_center(i+1), M_trans(i, i+1)];
    fillTri(x3, y3, M3);
end

% ---- D) Wall-edge quads: between last transition column and wall ----
% C-(i) hits wall at x_wall_pts(i+1). 
% The cell between adjacent C- waves near the wall:
%   corners: trans(i,N) [or center(i) if i==N], wall(i+1), wall(i+2), trans(i+1,N)[or center(i+1)]
% Actually: column of nodes just below wall = trans(:,N), except last C- (i=N) -> center(N)

% Get the "pre-wall" node for each C- wave i:
preWall_x = zeros(N,1);  preWall_y = zeros(N,1);  preWall_M = zeros(N,1);
for i = 1:N
    if i < N && ~isnan(x_trans(i,N))
        preWall_x(i) = x_trans(i,N);
        preWall_y(i) = y_trans(i,N);
        preWall_M(i) = M_trans(i,N);
    else
        preWall_x(i) = x_center(i);
        preWall_y(i) = y_center(i);
        preWall_M(i) = M_center(i);
    end
end

for i = 1:N-1
    % Quad: preWall(i), preWall(i+1), wall(i+2), wall(i+1)
    % wall indices: C-(i) -> wall_pts(i+1), C-(i+1) -> wall_pts(i+2)
    if i+2 > length(x_wall_pts); continue; end
    x4 = [preWall_x(i), preWall_x(i+1), x_wall_pts(i+2), x_wall_pts(i+1)];
    y4 = [preWall_y(i), preWall_y(i+1), y_wall_pts(i+2), y_wall_pts(i+1)];
    M4 = [preWall_M(i), preWall_M(i+1), M_wall(i+1),     M_wall(i)];
    fillQuad(x4, y4, M4);
end

% ---- E) First wall triangle: throat arc end -> first pre-wall node -> wall(1) ----
% The region between the throat arc endpoint, the first C- and the first wall point
x3 = [x_wall_pts(1), preWall_x(1), x_wall_pts(2)];
y3 = [y_wall_pts(1), preWall_y(1), y_wall_pts(2)];
M3 = [meyerMach(nu0(N),gamma), preWall_M(1), M_wall(1)];
fillTri(x3, y3, M3);

% ---- F) Exit uniform region (if wall was extended beyond last MOC point) ----
if length(x_wall_pts) > N+1
    xW1 = x_wall_pts(N+1); yW1 = y_wall_pts(N+1);
    xW2 = x_wall_pts(end); yW2 = y_wall_pts(end);
    x4 = [xW1, xW2, xW2, xW1];
    y4 = [yW1, yW2, 0,   0  ];
    patch(x4, y4, M_e, 'EdgeColor','none','FaceColor', getColor(M_e));
end

% -- Region 3: Wall strips (last transition row to wall points) --
for i = 1:N
    % Wall point index is i+1 in x_wall_pts
    xW = x_wall_pts(i+1); yW = y_wall_pts(i+1); MW = M_wall(i);
    xWprev = x_wall_pts(i); yWprev = y_wall_pts(i);

    if i > 1
        MWprev = M_wall(i-1);
    else
        MWprev = meyerMach(nu0(N), gamma);
    end
    
    % Previous interior node on this C-:
    if i < N && ~isnan(x_trans(i,N))
        xP = x_trans(i,N); yP = y_trans(i,N); MP = M_trans(i,N);
    elseif i == N
        xP = x_center(N); yP = y_center(N); MP = M_center(N);
    else
        xP = x_center(i); yP = y_center(i); MP = M_center(i);
    end
    
    % Triangle: wall(i), wall(i+1), interior node
    xs = [xWprev, xW, xP]; ys = [yWprev, yW, yP];
    Ms = [MWprev, MW, MP];
    M_p = mean(Ms);
    patch(xs, ys, M_p, 'EdgeColor','none','FaceColor',getColor(M_p));
end

% -- Fill the uniform exit region (last wall point to exit, if extended) --
if length(x_wall_pts) > N+1
    % The exit extension strip
    xW1 = x_wall_pts(N+1); yW1 = y_wall_pts(N+1);
    xW2 = x_wall_pts(end); yW2 = y_wall_pts(end);
    xs = [xW1, xW2, x_wall_pts(end), xW1];
    ys = [yW1, yW2, 0, 0];
    patch(xs, ys, M_e, 'EdgeColor','none','FaceColor',getColor(M_e));
end

% -------------------------------------------------------------------------
% Plot characteristic POLYLINES (connected, correct order, no zigzag)
% -------------------------------------------------------------------------
wave_color = [1 1 1];  % white lines on colored background, like reference image

% C+ waves: throat source -> centerline (downward kernel lines)
for j = 1:N
    poly = Cplus_poly{j};
    plot(poly(:,1), poly(:,2), '-', 'Color', [wave_color 0.7], 'LineWidth', 0.6);
end

% C- waves: centerline -> transition nodes -> wall (upward reflected lines)
for i = 1:N
    poly = Cminus_poly{i};
    plot(poly(:,1), poly(:,2), '-', 'Color', [wave_color 0.7], 'LineWidth', 0.6);
end

% -------------------------------------------------------------------------
% Wall and throat arc (drawn on top)
% -------------------------------------------------------------------------
delta_fine = linspace(0, delta_max, 100);
[xArc_fine, yArc_fine] = arcpoints(delta_fine);
plot(xArc_fine, yArc_fine, '-', 'Color', 'k', 'LineWidth', 2.5);

if length(x_wall_pts) > 2
    t_param = linspace(0, 1, length(x_wall_pts));
    t_fine  = linspace(0, 1, 200);
    x_wall_smooth = interp1(t_param, x_wall_pts, t_fine, 'pchip');
    y_wall_smooth = interp1(t_param, y_wall_pts, t_fine, 'pchip');
    plot(x_wall_smooth, y_wall_smooth, '-', 'Color', 'k', 'LineWidth', 2.5);
else
    plot(x_wall_pts, y_wall_pts, '-', 'Color', 'k', 'LineWidth', 2.5);
end

% Centerline
x_max_plot = max([all_x, x_wall_pts']) * 1.1;
plot([0 x_max_plot], [0 0], 'k--', 'LineWidth', 1);

% Exit plane
x_exit_line = x_wall_pts(end);
y_exit_line = y_wall_pts(end);
plot([x_exit_line x_exit_line], [0 y_exit_line], 'w:', 'LineWidth', 1.5);
text(x_exit_line*0.98, y_exit_line*0.5, sprintf('$M_e = %.1f$', M_e), ...
    'HorizontalAlignment', 'right', 'FontSize', 11, 'Color', 'black', 'Interpreter','latex');

% Node markers
plot(x0, y0, 'wo', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerSize', 4);
plot(x_center, y_center, 'w^', 'MarkerFaceColor', gray, 'MarkerSize', 4);
plot(x_wall_pts, y_wall_pts, 'ws', 'MarkerFaceColor', gray, 'MarkerSize', 5);

% -------------------------------------------------------------------------
% Colorbar, labels
% -------------------------------------------------------------------------
cb = colorbar;
cb.Label.String = '$M$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;
cb.TickLabelInterpreter = 'latex';
caxis([M_min_plot M_max_plot]);

xlabel('$x$ (m)', 'FontSize', 13);
ylabel('$y$ (m)', 'FontSize', 13);
title(sprintf('MOC Minimum Length Nozzle ($M_e = %.1f$, $N = %d$ waves, $\\gamma = %.2f$)', ...
    M_e, N, gamma), 'FontSize', 14);

hold off;

% -------------------------------------------------------------------------
% BONUS: Also plot the full nozzle (mirror about centerline) in a new figure
% -------------------------------------------------------------------------
plotFullNozzleProfile = false;
if(plotFullNozzleProfile)

figure('Name', 'Full Nozzle Profile', 'Position', [100 100 1400 500]);
hold on;
grid on;
axis equal;

% Upper wall (arc + computed wall)
plot(xArc_fine, yArc_fine, '-', 'Color', arizonaBlue, 'LineWidth', 2.5);
if length(x_wall_pts) > 2
    plot(x_wall_smooth, y_wall_smooth, '-', 'Color', arizonaRed, 'LineWidth', 2.5);
else
    plot(x_wall_pts, y_wall_pts, '-', 'Color', arizonaRed, 'LineWidth', 2.5);
end

% Lower wall (mirror)
plot(xArc_fine, -yArc_fine, '-', 'Color', arizonaBlue, 'LineWidth', 2.5);
if length(x_wall_pts) > 2
    plot(x_wall_smooth, -y_wall_smooth, '-', 'Color', arizonaRed, 'LineWidth', 2.5);
else
    plot(x_wall_pts, -y_wall_pts, '-', 'Color', arizonaRed, 'LineWidth', 2.5);
end

% Throat line
plot([0 0], [-y_wall_start y_wall_start], 'k-', 'LineWidth', 1.5);

% Exit plane
plot([x_exit_line x_exit_line], [-y_exit_line y_exit_line], ':', 'Color', river, 'LineWidth', 2);

% Centerline
plot([min(xArc_fine) x_exit_line*1.05], [0 0], 'k--', 'LineWidth', 0.5);

% Plot characteristics (upper half only) with Mach color gradient
for k = 1:size(char_lines_Cplus, 1)
    col = getColor(char_lines_Cplus(k, 5));
    plot([char_lines_Cplus(k,1) char_lines_Cplus(k,3)], ...
         [char_lines_Cplus(k,2) char_lines_Cplus(k,4)], ...
         '-', 'Color', [col 0.6], 'LineWidth', 0.8);
end
for k = 1:size(char_lines_Cminus, 1)
    col = getColor(char_lines_Cminus(k, 5));
    plot([char_lines_Cminus(k,1) char_lines_Cminus(k,3)], ...
         [char_lines_Cminus(k,2) char_lines_Cminus(k,4)], ...
         '--', 'Color', [col 0.6], 'LineWidth', 0.8);
end

colormap(cmap_custom);
cb2 = colorbar;
cb2.Label.String = 'Mach Number';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontSize = 13;
cb2.TickLabelInterpreter = 'latex';
caxis([M_min_plot M_max_plot]);

xlabel('$x$ (m)', 'FontSize', 13);
ylabel('$y$ (m)', 'FontSize', 13);
title(sprintf('Minimum Length Nozzle - Full Profile ($M_e = %.1f$)', M_e), 'FontSize', 14);

% Annotations
text(0, y_wall_start*1.2, sprintf('$A^*/2 = %.2f$ mm', h_t*1000), ...
    'HorizontalAlignment', 'center', 'FontSize', 10);
text(x_exit_line, y_exit_line*1.15, sprintf('$A_e/2 = %.2f$ mm', y_exit_line*1000), ...
    'HorizontalAlignment', 'center', 'FontSize', 10);

hold off;

end

fprintf('\n--- Plotting Complete ---\n');
fprintf('Wall contour points (x, y, delta):\n');
fprintf('%10s %10s %10s\n', 'x (mm)', 'y (mm)', 'delta (deg)');
for i = 1:length(x_wall_pts)
    fprintf('%10.3f %10.3f %10.2f\n', x_wall_pts(i)*1000, y_wall_pts(i)*1000, delta_wall_pts(i));
end