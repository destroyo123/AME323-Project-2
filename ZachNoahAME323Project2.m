clc; close all; clear;
%% AME 323 Project 2 - Supersonic Nozzle
%

%% Variables
G = 1.4;     % Gamma (air)
Me = 6.0;    % Design exit Mach number
N = 40;      % Number of characteristics (to be determined by convergence)


%% TASK 2 - PART A: Grid Convergence Study (Area Ratio Error Metric)

% Compute the EXACT isentropic area ratio for Mach 6
A_ratio_exact = (1/Me) * ((2/(G+1)) * (1 + (G-1)/2 * Me^2))^((G+1)/(2*(G-1)));

% Range of N values to test
N_values = 5:1:104;  % N from 5 to 104
error_AR = zeros(size(N_values));  % Area ratio error storage

cutoff = 0.001;       % Convergence threshold for area ratio error
N_optimal = N_values(end);  % Default if threshold never met
found_optimal = false;


for idx = 1:length(N_values)
    n = N_values(idx);
    
    % --- Run MOC for this N ---
    Cm = zeros(n,n);
    Cp = zeros(n,n);
    Theta = zeros(n,n);
    Mu = zeros(n,n);
    M = zeros(n,n);
    Nu = zeros(n,n);
    x = zeros(n,n);
    y = zeros(n,n);
    
    [~, NuMax, ~] = PMF(G, Me, 0, 0);
    ThetaMax = NuMax/2;
    dT = ThetaMax/n;
    ThetaArc = (0:dT:ThetaMax)';
    NuArc = ThetaArc;
    CmArc = ThetaArc + NuArc;
    MuArc = zeros(size(ThetaArc));
    
    for i = 1:length(ThetaArc)
        [~, ~, MuArc(i)] = PMF(G, 0, NuArc(i), 0);
    end
    
    y0 = 1;
    ThroatCurveRadius = 1.5*y0;
    [xarc, yarc] = Arc(ThroatCurveRadius, ThetaArc);
    yarc = yarc + y0;
    
    Cm(:,1) = CmArc(2:n+1);
    Theta(:,1) = ThetaArc(2:n+1);
    Nu(:,1) = Theta(:,1);
    Cp(:,1) = Theta(:,1) - Nu(:,1);
    M(1,1) = 1.0001;
    Nu(1,1) = 0;
    Mu(1,1) = 90;
    y(1,1) = 0;
    x(1,1) = xarc(2) + (y(1,1)-yarc(2))/tand((ThetaArc(2)-MuArc(2)+ThetaArc(2)-MuArc(2))/2);
    
    for i = 2:n
        [M(i,1), Nu(i,1), Mu(i,1)] = PMF(G, 0, Nu(i,1), 0);
        s1 = tand((ThetaArc(i+1) - MuArc(i+1) + Theta(i,1) - Mu(i,1))/2);
        s2 = tand((Theta(i-1,1) + Mu(i-1,1) + Theta(i,1) + Mu(i,1))/2);
        x(i,1) = ((y(i-1,1) - x(i-1,1)*s2) - (yarc(i+1) - xarc(i+1)*s1))/(s1-s2);
        y(i,1) = y(i-1,1) + (x(i,1) - x(i-1,1))*s2;
    end
    
    for j = 2:n
        for i = 1:n+1-j
            Cm(i,j) = Cm(i+1,j-1);
            if i == 1
                Theta(i,j) = 0;
                Cp(i,j) = -Cm(i,j);
                Nu(i,j) = Cm(i,j);
                [M(i,j), Nu(i,j), Mu(i,j)] = PMF(G, 0, Nu(i,j), 0);
                s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
                x(i,j) = x(i+1,j-1) - y(i+1,j-1)/s1;
                y(i,j) = 0;
            else
                Cp(i,j) = Cp(i-1,j);
                Theta(i,j) = (Cm(i,j)+Cp(i,j))/2;
                Nu(i,j) = (Cm(i,j)-Cp(i,j))/2;
                [M(i,j), Nu(i,j), Mu(i,j)] = PMF(G, 0, Nu(i,j), 0);
                s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
                s2 = tand((Theta(i-1,j)+Mu(i-1,j)+Theta(i,j)+Mu(i,j))/2);
                x(i,j) = ((y(i-1,j)-x(i-1,j)*s2)-(y(i+1,j-1)-x(i+1,j-1)*s1))/(s1-s2);
                y(i,j) = y(i-1,j) + (x(i,j)-x(i-1,j))*s2;
            end
        end
    end
    
    % --- Compute wall contour to get exit area ratio ---
    xwall = zeros(2*n,1); ywall = xwall; ThetaWall = ywall;
    xwall(1:n) = xarc(2:n+1); ywall(1:n) = yarc(2:n+1);
    ThetaWall(1:n) = ThetaArc(2:n+1);
    
    for i = 1:n-1
        ThetaWall(n+i) = ThetaWall(n-i);
    end
    
    for i = 1:n
        s1 = tand((ThetaWall(n+i-1)+ThetaWall(n+i))/2);
        s2 = tand(Theta(n+1-i,i)+Mu(n+1-i,i));
        xwall(n+i) = ((y(n+1-i,i)-x(n+1-i,i)*s2)-(ywall(n+i-1)-xwall(n+i-1)*s1))/(s1-s2);
        ywall(n+i) = ywall(n+i-1) + (xwall(n+i)-xwall(n+i-1))*s1;
    end
    
    % --- EXIT AREA RATIO (direct geometric output) ---
    y_exit = ywall(end);
    A_ratio_computed = y_exit / y0;  % A_exit / A* for 2D planar nozzle
    
    % --- AREA RATIO ERROR (primary convergence metric) ---
    error_AR(idx) = abs(A_ratio_computed - A_ratio_exact) / A_ratio_exact;
    
    % --- Check convergence ---
    if error_AR(idx) < cutoff && ~found_optimal
        N_optimal = n;
        found_optimal = true;
    end
end

fprintf("===== TASK 2 - PART A: Grid Convergence Study =====\n")
fprintf(" Exact isentropic area ratio for M = %.1f: A/A* = %.6f\n", Me, A_ratio_exact)
fprintf("\n===== TASK 2 - PART B: Justification of Number of Lines =====\n")
fprintf(" Optimized N = %d\n", N_optimal)
fprintf(" This was determined by requiring the relative error in\n")
fprintf(" exit area ratio to be less than %.2e\n", cutoff)
if found_optimal
    opt_idx = find(N_values == N_optimal, 1);
    fprintf(" At N = %d: A/A* computed = %.6f, Relative Error = %.2e\n", ...
        N_optimal, ...
        (1/Me)*((2/(G+1))*(1+(G-1)/2*Me^2))^((G+1)/(2*(G-1))) * (1 + error_AR(opt_idx)), ...
        error_AR(opt_idx))
end
%% PART A PLOT - Area Ratio Error vs N
figure('Position', [0 700 500 250]);

% Main plot: relative error in area ratio
plot(N_values, error_AR, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'r');
hold on;

% Mark optimal N
if found_optimal
    opt_idx = find(N_values == N_optimal, 1);
    semilogy(N_optimal, error_AR(opt_idx), 'bp', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end

set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel('Number of Characteristic Lines, N', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Relative Error in Exit Area Ratio', ...
    'FontName', 'Times New Roman', 'FontSize', 13);
title('TASK 2 - Part A: Grid Convergence Study (Area Ratio Error Metric)', ...
    'FontName', 'Times New Roman', 'FontSize', 13);
grid on;

legend({'Relative Area Ratio Error', sprintf('Optimal N = %d', N_optimal)}, 'FontName', 'Times New Roman', 'Location', 'northeast');

hold off;

%% TASK 2 - PARTS B and C
% Initialize Variables
Cm = zeros(N,N); 
Cp = zeros(N,N); 
Theta = zeros(N,N);
Mu = zeros(N,N); 
M = zeros(N,N); 
Nu = zeros(N,N);
x = zeros(N,N); 
y = zeros(N,N);

Nu_center = [];
p_p0center = [];
T_T0center = [];
M_center = [];
x_center = [];

% Nu for Mach 6
[~, NuMax, ~] = PMF(G, Me, 0, 0);
ThetaMax = NuMax/2; dT = ThetaMax/N;


ThetaArc = (0:dT:ThetaMax)';
NuArc = ThetaArc; 
CmArc = ThetaArc + NuArc;
MuArc = zeros(size(ThetaArc));

for i = 1:length(ThetaArc)
    [~, ~, MuArc(i)] = PMF(G, 0, NuArc(i), 0);
end

y0 = 1; % For y/y0, at the throat this should be 1
ThroatCurveRadius = 1.5*y0;
[xarc, yarc] = Arc(ThroatCurveRadius, ThetaArc);
yarc = yarc + y0;

Cm(:,1) = CmArc(2:N+1); 
Theta(:,1) = ThetaArc(2:N+1);
Nu(:,1) = Theta(:,1); 
Cp(:,1) = Theta(:,1) - Nu(:,1);
M(1,1) = 1.0001; Nu(1,1) = 0; 
Mu(1,1) = 90; 
y(1,1) = 0;
x(1,1) = xarc(2) + (y(1,1)-yarc(2))/tand((ThetaArc(2)-MuArc(2)+ThetaArc(2)-MuArc(2))/2);

for i=2:N
    [M(i,1), Nu(i,1), Mu(i,1)] = PMF(G, 0, Nu(i,1), 0);
    s1 = tand((ThetaArc(i+1)-MuArc(i+1)+Theta(i,1)-Mu(i,1))/2);
    s2 = tand((Theta(i-1,1)+Mu(i-1,1)+Theta(i,1)+Mu(i,1))/2);
    x(i,1) = ((y(i-1,1)-x(i-1,1)*s2)-(yarc(i+1)-xarc(i+1)*s1))/(s1-s2);
    y(i,1) = y(i-1,1)+(x(i,1)-x(i-1,1))*s2;
end

for j=2:N
    for i=1:N+1-j
        Cm(i,j) = Cm(i+1,j-1);
        if i==1
            Theta(i,j)=0; 
            Cp(i,j)=-Cm(i,j); 
            Nu(i,j)=Cm(i,j);
            [M(i,j),Nu(i,j),Mu(i,j)] = PMF(G,0,Nu(i,j),0);
            s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
            x(i,j) = x(i+1,j-1)-y(i+1,j-1)/s1; 
            y(i,j)=0;
            
            Nu_center(end + 1) = Nu(i,j);
            
            fun = @(Mach) rad2deg(sqrt((G+1)/(G-1)) * atan(sqrt((G-1)/(G+1)*(Mach^2-1))) - atan(sqrt(Mach^2-1))) - Nu_center(end);
            M_center(end + 1) = fzero(fun, [1.0001 30]);

            p_p0center(end + 1) = (1 + ((G - 1)/2)*(M_center(end))^2)^(-G/(G-1));
            T_T0center(end + 1) = (1 + ((G - 1)/2)*(M_center(end))^2)^(-1);
            x_center(end + 1) = x(i,j);

        else
            Cp(i,j)=Cp(i-1,j);
            Theta(i,j)=(Cm(i,j)+Cp(i,j))/2; Nu(i,j)=(Cm(i,j)-Cp(i,j))/2;
            [M(i,j),Nu(i,j),Mu(i,j)] = PMF(G,0,Nu(i,j),0);
            s1 = tand((Theta(i+1,j-1)-Mu(i+1,j-1)+Theta(i,j)-Mu(i,j))/2);
            s2 = tand((Theta(i-1,j)+Mu(i-1,j)+Theta(i,j)+Mu(i,j))/2);
            x(i,j) = ((y(i-1,j)-x(i-1,j)*s2)-(y(i+1,j-1)-x(i+1,j-1)*s1))/(s1-s2);
            y(i,j) = y(i-1,j)+(x(i,j)-x(i-1,j))*s2;
        end
    end
end

xwall = zeros(2*N,1); 
ywall = xwall; 
ThetaWall = ywall;
xwall(1:N) = xarc(2:N+1); 
ywall(1:N) = yarc(2:N+1);
ThetaWall(1:N) = ThetaArc(2:N+1);

for i=1:N-1
    ThetaWall(N+i)=ThetaWall(N-i);
end

for i=1:N
    s1=tand((ThetaWall(N+i-1)+ThetaWall(N+i))/2);
    s2=tand(Theta(N+1-i,i)+Mu(N+1-i,i));
    xwall(N+i)=((y(N+1-i,i)-x(N+1-i,i)*s2)-(ywall(N+i-1)-xwall(N+i-1)*s1))/(s1-s2);
    ywall(N+i)=ywall(N+i-1)+(xwall(N+i)-xwall(N+i-1))*s1;
end

%% PARTS B and C PLOTS
% Plotting Nozzle Contour w/ Characteristic Lines, Plotting p/p0, T/T0, and M

figure('Position', [0 400 500 250]);
    
    plot(xwall,ywall,'-','LineWidth',2); axis equal; hold on
    plot(xarc,yarc,'k-','LineWidth',1.5);

    for i=1:N-1
        plot(x(1:N+1-i,i), y(1:N+1-i,i)); 
    end

    for i=1:N
        plot([xarc(i) x(i,1)],[yarc(i) y(i,1)]);
        plot([x(N+1-i,i) xwall(i+N)],[y(N+1-i,i) ywall(i+N)]);
    end

    for c=1:N
        for r=2:N+1-c
            plot([x(c,r) x(c+1,r-1)],[y(c,r) y(c+1,r-1)]);
        end 
    end

    xlabel('Length [x/y0]', 'FontName', "Times New Roman"); 
    ylabel('Height [y/y0]', 'FontName', "Times New Roman");
    title('=== TASK 2 - PART B: Mach 6 Nozzle (Optimized N) ===', 'FontName', "Times New Roman"); 
    axes = gca;
    set(axes, 'FontName', 'Times New Roman');

    grid on

figure('Position', [0 100 500 250]);
    plot(p_p0center, M_center, 'LineWidth', 2);
    hold on;
    plot(T_T0center, M_center, 'LineWidth', 2);
    xlabel('Mach Number', 'FontName', "Times New Roman"); 
    ylabel('Ratios', 'FontName', "Times New Roman");
    title('=== TASK 2 - PART C: P/P0 and T/T0 vs Mach Number Along Centerline ===', 'FontName', "Times New Roman"); 
    axes = gca;
    set(axes, 'FontName', 'Times New Roman');
    legend({'P/P0', 'T/T0'}, 'FontName', 'Times New Roman');
    xlim([-0.05 0.85]);

    hold off
     
P_cr3 = (1 + (G-1)/2 * Me^2)^(-G/(G-1));

fprintf('\n===== TASK 2 - PART C: Critical Pressure Ratio =====\n');
fprintf(' P_cr3 (p_back/p_0) for shock-free operation at Mach %.1f: %.6e\n', Me, P_cr3);

%% TASK 2 - PART D
% Collect scattered data points inside nozzle
X_all = [];
Y_all = [];
Mach_all = [];
P_ratio_all = [];
T_ratio_all = [];

% 1) Internal MOC grid points
for j = 1:N
    for i = 1:N+1-j
        if x(i,j) ~= 0 || y(i,j) ~= 0 || (i==1 && j==1)
            X_all(end+1) = x(i,j);
            Y_all(end+1) = y(i,j);
            Mach_all(end+1) = M(i,j);
            P_ratio_all(end+1) = (1 + (G-1)/2*M(i,j)^2)^(-G/(G-1));
            T_ratio_all(end+1) = (1 + (G-1)/2*M(i,j)^2)^(-1);
        end
    end
end

% 2) Arc (expansion) wall points
for i = 2:N+1
    [M_arc, ~, ~] = PMF(G, 0, NuArc(i), 0);
    X_all(end+1) = xarc(i);
    Y_all(end+1) = yarc(i);
    Mach_all(end+1) = M_arc;
    P_ratio_all(end+1) = (1 + (G-1)/2*M_arc^2)^(-G/(G-1));
    T_ratio_all(end+1) = (1 + (G-1)/2*M_arc^2)^(-1);
end

% 3) Straightening section wall points
for i = 1:N
    idx_r = N+1-i;
    M_wall_pt = M(idx_r, i);
    X_all(end+1) = xwall(N+i);
    Y_all(end+1) = ywall(N+i);
    Mach_all(end+1) = M_wall_pt;
    P_ratio_all(end+1) = (1 + (G-1)/2*M_wall_pt^2)^(-G/(G-1));
    T_ratio_all(end+1) = (1 + (G-1)/2*M_wall_pt^2)^(-1);
end

% 4) Throat points
X_all(end+1) = 0; Y_all(end+1) = y0;
Mach_all(end+1) = 1.0;
P_ratio_all(end+1) = (1 + (G-1)/2*1.0^2)^(-G/(G-1));
T_ratio_all(end+1) = (1 + (G-1)/2*1.0^2)^(-1);

X_all(end+1) = 0; Y_all(end+1) = 0;
Mach_all(end+1) = 1.0;
P_ratio_all(end+1) = (1 + (G-1)/2*1.0^2)^(-G/(G-1));
T_ratio_all(end+1) = (1 + (G-1)/2*1.0^2)^(-1);

% ============================================================
% 5) ADD UNIFORM FLOW REGION BEYOND LAST EXPANSION WAVE
% ============================================================
% After the last characteristic, flow is uniform at Me
M_uniform = Me;
P_uniform = (1 + (G-1)/2*M_uniform^2)^(-G/(G-1));
T_uniform = (1 + (G-1)/2*M_uniform^2)^(-1);

% Determine the x-extent of the uniform region
x_last_char = max(X_all);       % rightmost MOC point
x_exit = xwall(end);            % nozzle exit
y_exit = ywall(end);

% Extend slightly beyond the exit for coverage
x_extend = x_exit + 0.1*(x_exit - x_last_char);

% Create a grid of uniform-flow points from the last characteristic to exit
n_uniform_x = 20;
n_uniform_y = 10;
x_uni_vec = linspace(x_last_char, x_extend, n_uniform_x);
y_uni_vec = linspace(0, y_exit, n_uniform_y);

for ii = 1:n_uniform_x
    for jj = 1:n_uniform_y
        % Interpolate wall height at this x to stay inside nozzle
        y_wall_here = interp1([xwall(N); xwall(N+1:2*N)], ...
                              [ywall(N); ywall(N+1:2*N)], ...
                              x_uni_vec(ii), 'linear', y_exit);
        if y_uni_vec(jj) <= y_wall_here
            X_all(end+1) = x_uni_vec(ii);
            Y_all(end+1) = y_uni_vec(jj);
            Mach_all(end+1) = M_uniform;
            P_ratio_all(end+1) = P_uniform;
            T_ratio_all(end+1) = T_uniform;
        end
    end
end

% ============================================================
% Build full wall boundary for masking
% ============================================================
xwall_full = [0; xarc(2:N+1); xwall(N+1:2*N)];
ywall_full = [y0; yarc(2:N+1); ywall(N+1:2*N)];
[xwall_sorted, sort_idx] = sort(xwall_full);
ywall_sorted = ywall_full(sort_idx);
[xwall_unique, unique_idx] = unique(xwall_sorted);
ywall_unique = ywall_sorted(unique_idx);

% Create interpolation grid
xq = linspace(min(X_all), max(X_all), 400);
yq = linspace(0, max(Y_all)*1.05, 200);
[Xq, Yq] = meshgrid(xq, yq);

% Remove duplicate points before interpolation
coords = [X_all(:), Y_all(:)];
[unique_coords, ia, ~] = unique(coords, 'rows', 'stable');
X_unique = unique_coords(:,1);
Y_unique = unique_coords(:,2);
Mach_unique = Mach_all(ia)';
P_unique = P_ratio_all(ia)';
T_unique = T_ratio_all(ia)';

% Create interpolants without duplicate warnings
F_M = scatteredInterpolant(X_unique, Y_unique, Mach_unique, 'natural', 'none');
F_P = scatteredInterpolant(X_unique, Y_unique, P_unique, 'natural', 'none');
F_T = scatteredInterpolant(X_unique, Y_unique, T_unique, 'natural', 'none');

Mach_q = F_M(Xq, Yq);
P_q    = F_P(Xq, Yq);
T_q    = F_T(Xq, Yq);

% Mask points outside nozzle boundary
ywall_interp = interp1(xwall_unique, ywall_unique, xq, 'linear', NaN);
for k = 1:length(yq)
    for m = 1:length(xq)
        if isnan(ywall_interp(m)) || Yq(k,m) > ywall_interp(m) || Xq(k,m) < 0
            Mach_q(k,m) = NaN;
            P_q(k,m)    = NaN;
            T_q(k,m)    = NaN;
        end
    end
end

% ============================================================
% SUBPLOT LAYOUT — all three contours in one figure
% ============================================================
figure('Position', [50 50 900 900]);

% --- Subplot 1: Mach Number ---
ax1 = subplot(3,1,1);
contourf(ax1, Xq, Yq, Mach_q, 30, 'LineColor', 'none');
hold(ax1, 'on');
plot(ax1, xwall_full, ywall_full, 'k', 'LineWidth', 2);
plot(ax1, xarc, yarc, 'k', 'LineWidth', 2);
colormap(ax1, jet);          % colormap applied to THIS axes only
cb1 = colorbar(ax1);
ylabel(cb1, 'Mach Number', 'FontName', 'Times New Roman');
xlabel(ax1, 'x / y_0', 'FontName', 'Times New Roman');
ylabel(ax1, 'y / y_0', 'FontName', 'Times New Roman');
title(ax1, 'TASK 2 - PART D: Mach Number Contours', 'FontName', 'Times New Roman');
set(ax1, 'FontName', 'Times New Roman');
axis(ax1, 'equal', 'tight');
grid(ax1, 'on');

% --- Subplot 2: Pressure Ratio ---
ax2 = subplot(3,1,2);
contourf(ax2, Xq, Yq, P_q, 30, 'LineColor', 'none');
hold(ax2, 'on');
plot(ax2, xwall_full, ywall_full, 'k', 'LineWidth', 2);
plot(ax2, xarc, yarc, 'k', 'LineWidth', 2);
colormap(ax2, jet);          % separate colormap for this axes
cb2 = colorbar(ax2);
ylabel(cb2, 'p / p_0', 'FontName', 'Times New Roman');
xlabel(ax2, 'x / y_0', 'FontName', 'Times New Roman');
ylabel(ax2, 'y / y_0', 'FontName', 'Times New Roman');
title(ax2, 'TASK 2 - PART D: Pressure Ratio (p/p_0) Contours', 'FontName', 'Times New Roman');
set(ax2, 'FontName', 'Times New Roman');
axis(ax2, 'equal', 'tight');
grid(ax2, 'on');

% --- Subplot 3: Temperature Ratio ---
ax3 = subplot(3,1,3);
contourf(ax3, Xq, Yq, T_q, 30, 'LineColor', 'none');
hold(ax3, 'on');
plot(ax3, xwall_full, ywall_full, 'k', 'LineWidth', 2);
plot(ax3, xarc, yarc, 'k', 'LineWidth', 2);
colormap(ax3, jet);          % separate colormap for this axes
cb3 = colorbar(ax3);
ylabel(cb3, 'T / T_0', 'FontName', 'Times New Roman');
xlabel(ax3, 'x / y_0', 'FontName', 'Times New Roman');
ylabel(ax3, 'y / y_0', 'FontName', 'Times New Roman');
title(ax3, 'TASK 2 - PART D: Temperature Ratio (T/T_0) Contours', 'FontName', 'Times New Roman');
set(ax3, 'FontName', 'Times New Roman');
axis(ax3, 'equal', 'tight');
grid(ax3, 'on');

%% Task 2 - Part E
% Physical scaling based on 350 mm exit height
h_exit = 0.350; % Exit height [m] (350 mm)
y0_physical = (h_exit/2) / ywall(end); % Physical throat half-height [m]
L_nozzle_physical = xwall(end) * y0_physical; % Physical nozzle length [m]

% Uniform flow begins at the last centerline characteristic intersection
x_uniform_start_physical = x_center(end) * y0_physical; % [m]
UniformFlowLength = (xwall(end) - x_center(end)) * y0_physical; % [m]

fprintf("\n\n==== TASK 2 - PART E: Uniform Flow Start & Length ====\n")
fprintf(" y0 (physical throat half-height) = %.4f mm\n", y0_physical*1000)
fprintf(" Total nozzle expansion length = %.4f m = %.2f mm\n", L_nozzle_physical, L_nozzle_physical*1000)
fprintf(" Uniform Flow at Mach 6 begins at x = %.4f m = %.2f mm\n", x_uniform_start_physical, x_uniform_start_physical*1000)
fprintf(" Uniform flow region length = %.4f m = %.2f mm\n", UniformFlowLength, UniformFlowLength*1000)



%% TASK 3: Estimate Tunnel Characteristics

%% Given Parameters
G = 1.4;            % Gamma (air)
Me = 6.0;           % Exit Mach number
R_air = 287;        % Specific gas constant for air [J/(kg·K)]
h_exit = 0.350;     % Exit height [m] (350 mm)
w = 0.500;          % Nozzle width [m] (500 mm)
p_b = 0.1;          % Back pressure [psia]
V_tank = 40;        % Tank volume [m^3]
T_tank = 300;       % Tank temperature [K] (constant due to thermal matrix)

% Sutherland's Law constants
mu_ref = 1.716e-5;  % [Pa·s]
T_ref = 273.15;     % [K]
S = 110.4;          % [K]

%% P_cr3 - Critical Pressure Ratio for Shock-Free Operation
P_cr3 = (1 + (G-1)/2 * Me^2)^(-G/(G-1));
fprintf('\n\n===== TASK 3 =====\n');
fprintf('P_cr3 (p_exit/p_0) for shock-free operation at Mach %.1f: %.6e\n', Me, P_cr3);

%% Minimum Total Temperature (to avoid liquefaction)
% T_static >= 55 K at freestream
% T/T0 = (1 + (G-1)/2 * M^2)^(-1)
T_min_static = 55; % [K]
T_T0_exit = (1 + (G-1)/2 * Me^2)^(-1);
T0_min = T_min_static / T_T0_exit;

fprintf('\nMinimum Total Temperature:\n');
fprintf('  T/T0 at Mach %.1f = %.6f\n', Me, T_T0_exit);
fprintf('  T0_min = %.2f K to keep T_static >= 55 K\n', T0_min);

%% Minimum Total Pressure for Steady Operation (no shocks at exit)
% For no shocks at the exit, the back pressure must equal the exit static pressure:
%   p_b = p_exit = p_0 * P_cr3
% Therefore: p_0_min = p_b / P_cr3

p_b_Pa = p_b * 6894.76;  % Convert psia to Pa
p0_min_Pa = p_b_Pa / P_cr3;
p0_min_psia = p0_min_Pa / 6894.76;

fprintf('\nMinimum Total Pressure for Steady Operation (no shocks at exit):\n');
fprintf('  p_0_min = p_b / P_cr3 = %.2f / %.6e = %.2f psia\n', p_b, P_cr3, p0_min_psia);
fprintf('  p_0_min = %.2f Pa = %.4f kPa\n', p0_min_Pa, p0_min_Pa/1000);

%% Mass Flow Rate vs. Total Pressure
% A* = throat area. For 2D planar nozzle: A_exit/A* = area ratio
% A_exit = h_exit * w
% A* = A_exit / A_ratio_exact

A_ratio_exact = (1/Me) * ((2/(G+1)) * (1 + (G-1)/2 * Me^2))^((G+1)/(2*(G-1)));
A_exit = h_exit * w;          % [m^2]
A_star = A_exit / A_ratio_exact; % Throat area [m^2]

fprintf('\nGeometric Parameters:\n');
fprintf('  A/A* (exact) = %.4f\n', A_ratio_exact);
fprintf('  A_exit = %.4f m^2\n', A_exit);
fprintf('  A* (throat area) = %.6f m^2\n', A_star);

% Mass flow rate: mdot = p0 * A* * sqrt(G/(R*T0)) * (2/(G+1))^((G+1)/(2*(G-1)))
% Use T0 = T0_min for mass flow calculations

T0 = T0_min; % Use minimum total temperature

p0_range_psia = linspace(p0_min_psia, 500, 200);  % Operating range
p0_range_Pa = p0_range_psia * 6894.76;

% Isentropic mass flow rate through choked throat
mdot = p0_range_Pa * A_star .* sqrt(G / (R_air * T0)) * ((2/(G+1))^((G+1)/(2*(G-1))));

figure('Position', [50 550 560 420]);
plot(p0_range_psia, mdot, 'b-', 'LineWidth', 2);
xlabel('Total Pressure, p_0 [psia]', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Mass Flow Rate [kg/s]', 'FontName', 'Times New Roman', 'FontSize', 12);
title('TASK 3: Mass Flow Rate vs. Total Pressure', 'FontName', 'Times New Roman', 'FontSize', 13);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);

fprintf('\nMass Flow Rate Range:\n');
fprintf('  At p0 = %.2f psia (min): mdot = %.4f kg/s\n', p0_min_psia, mdot(1));
fprintf('  At p0 = 500 psia (max):  mdot = %.4f kg/s\n', 500, mdot(end));

%% Run Time vs. Total Pressure
% Tank: V = 40 m^3, T = 300 K (isothermal due to thermal matrix)
% Mass in tank: m_tank = p_tank * V / (R * T_tank)
% Tunnel runs until tank pressure drops to p0 (regulated pressure)
% Usable mass: m_usable = (p_tank_initial - p0) * V / (R * T_tank)
% But the regulator supplies CONSTANT p0 from the tank, so:
% Available mass = (p_max_tank - p0) * V / (R * T_tank)
% Actually, the tank starts at some pressure and the regulator holds p0 constant
% until tank pressure drops below p0.
% 
% Initial tank pressure: we need to determine this. The problem says air is 
% "stored at 300 K" in a 40 m^3 tank. The maximum regulated pressure is 500 psia.
% The tank must be at least at p0 to supply the regulator.
% 
% Assuming tank starts at max pressure (500 psia) and regulator holds p0 constant:
% Usable mass = (p_tank_max - p0) * V / (R * T_tank)
% Run time = usable mass / mdot

p_tank_max_Pa = 500 * 6894.76;  % Maximum tank pressure (500 psia)

run_time = zeros(size(p0_range_Pa));
for k = 1:length(p0_range_Pa)
    p0_k = p0_range_Pa(k);
    m_usable = (p_tank_max_Pa - p0_k) * V_tank / (R_air * T_tank);
    mdot_k = p0_k * A_star * sqrt(G / (R_air * T0)) * ((2/(G+1))^((G+1)/(2*(G-1))));
    run_time(k) = m_usable / mdot_k;
end

figure('Position', [650 550 560 420]);
plot(p0_range_psia, run_time, 'r-', 'LineWidth', 2);
xlabel('Total Pressure, p_0 [psia]', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Run Time [s]', 'FontName', 'Times New Roman', 'FontSize', 12);
title('TASK 3: Tunnel Run Time vs. Total Pressure', 'FontName', 'Times New Roman', 'FontSize', 13);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);

fprintf('\nRun Time:\n');
fprintf('  At p0 = %.2f psia (min): t_run = %.2f s\n', p0_min_psia, run_time(1));
fprintf('  At p0 = 500 psia (max):  t_run = %.2f s\n', 500, run_time(end));

%% Unit Reynolds Number vs. Total Pressure
% Re' = rho * U / mu [1/m]
% At freestream (exit): 
%   T_exit = T0 * T_T0_exit
%   p_exit = p0 * P_cr3
%   rho_exit = p_exit / (R * T_exit)
%   a_exit = sqrt(G * R * T_exit)
%   U_exit = Me * a_exit
%   mu from Sutherland's law at T_exit

T_exit = T0_min * T_T0_exit;  % Static temperature at exit [K]
a_exit = sqrt(G * R_air * T_exit);  % Speed of sound at exit [m/s]
U_exit = Me * a_exit;               % Flow velocity at exit [m/s]

% Sutherland's law for viscosity at T_exit
mu_exit = mu_ref * (T_exit/T_ref)^(3/2) * (T_ref + S) / (T_exit + S);

% Unit Reynolds number as function of p0
p_exit_Pa = p0_range_Pa * P_cr3;
rho_exit = p_exit_Pa / (R_air * T_exit);
Re_unit = rho_exit * U_exit / mu_exit;  % [1/m]

figure('Position', [650 100 560 420]);
plot(p0_range_psia, Re_unit/1e6, 'g-', 'LineWidth', 2);
xlabel('Total Pressure, p_0 [psia]', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Unit Reynolds Number [×10^6 /m]', 'FontName', 'Times New Roman', 'FontSize', 12);
title('TASK 3: Unit Reynolds Number vs. Total Pressure', 'FontName', 'Times New Roman', 'FontSize', 13);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);

fprintf('\nUnit Reynolds Number:\n');
fprintf('  T_exit = %.2f K\n', T_exit);
fprintf('  U_exit = %.2f m/s\n', U_exit);
fprintf('  mu_exit = %.4e Pa·s\n', mu_exit);
fprintf('  At p0 = %.2f psia: Re'' = %.4e /m\n', p0_min_psia, Re_unit(1));
fprintf('  At p0 = 500 psia:  Re'' = %.4e /m\n', 500, Re_unit(end));

%% Additional Geometric Characteristics

% Physical scaling: normalized exit half-height maps to h_exit/2 = 175 mm
y0_physical = (h_exit/2) / ywall(end); % Physical throat half-height [m]

% 1) Nozzle expansion length (from throat to exit)
L_nozzle_physical = xwall(end) * y0_physical; % Physical nozzle length [m]

% Throat dimensions
h_throat = 2 * y0_physical; % Full throat height [m]
A_star_check = h_throat * w; % Throat area [m^2]

fprintf('\nNozzle Physical Dimensions (scaled to h_exit = 350 mm):\n');
fprintf(' Normalized exit half-height (ywall(end)/y0) = %.4f\n', ywall(end));
fprintf(' Physical throat half-height (y0) = %.4f mm\n', y0_physical*1000);
fprintf(' Physical throat full height = %.4f mm\n', h_throat*1000);
fprintf(' Throat area (A*) = %.6f m^2 = %.2f mm^2\n', A_star_check, A_star_check*1e6);
fprintf(' Nozzle expansion length = %.4f m = %.2f mm\n', L_nozzle_physical, L_nozzle_physical*1000);
fprintf(' Exit height = %.1f mm (as specified)\n', h_exit*1000);
fprintf(' Exit area = %.4f m^2\n', A_exit);
fprintf(' Area ratio (A_exit/A*) = %.4f (exact = %.4f)\n', h_exit/h_throat, A_ratio_exact);

% 2) Settling chamber cross-sectional area for M <= 0.05
M_sc = 0.05;
A_ratio_sc = (1/M_sc) * ((2/(G+1)) * (1 + (G-1)/2 * M_sc^2))^((G+1)/(2*(G-1)));
A_sc = A_star_check * A_ratio_sc; % Settling chamber area [m^2]
h_sc = A_sc / w; % Settling chamber height [m]

fprintf('\nSettling Chamber:\n');
fprintf(' M_sc = %.2f\n', M_sc);
fprintf(' A_sc/A* = %.4f\n', A_ratio_sc);
fprintf(' A_sc = %.4f m^2\n', A_sc);
fprintf(' h_sc = %.4f m = %.2f mm\n', h_sc, h_sc*1000);

%% Nozzle Forces

% Thrust: F_thrust = mdot * U_exit + (p_exit - p_b) * A_exit
% At a given p0, compute thrust

% Use p0 = p0_min for a representative calculation (or you can do full range)
p0_force = p0_min_Pa;  % Use minimum operating pressure
mdot_force = p0_force * A_star * sqrt(G / (R_air * T0)) * ((2/(G+1))^((G+1)/(2*(G-1))));
p_exit_force = p0_force * P_cr3;

% At design condition, p_exit = p_b (matched), so pressure thrust term ≈ 0
% But let's compute it generally:
F_thrust = mdot_force * U_exit + (p_exit_force - p_b_Pa) * A_exit;

fprintf('\nNozzle Forces (at p0_min = %.2f psia):\n', p0_min_psia);
fprintf('  mdot = %.4f kg/s\n', mdot_force);
fprintf('  Thrust (axial) = %.2f N\n', F_thrust);

% Vertical force on single nozzle plane (pressure integration on expansion wall)
% This requires integrating pressure along the wall contour
% p_wall(x) from the MOC solution

% Compute wall pressure distribution
p_wall = zeros(2*N, 1);

% Expansion section (arc region): Mach varies along arc
for i = 1:N
    [M_w, ~, ~] = PMF(G, 0, NuArc(i+1), 0);
    p_wall(i) = p0_force * (1 + (G-1)/2 * M_w^2)^(-G/(G-1));
end

% Straightening section: use Mach from last characteristic hitting wall
for i = 1:N
    idx_r = N+1-i;
    M_w = M(idx_r, i);
    p_wall(N+i) = p0_force * (1 + (G-1)/2 * M_w^2)^(-G/(G-1));
end

% Scale wall coordinates to physical dimensions
xwall_phys = xwall * y0_physical;  % [m]
ywall_phys = ywall * y0_physical;  % [m]

% Vertical force = integral of p * dx * w (pressure acting on horizontal projection)
% F_vertical = w * integral(p dA_y) where dA_y = dx * w for projection
% More precisely: F_y = w * sum(p * dx) for the vertical component

F_vertical = 0;
for i = 2:2*N
    dx = xwall_phys(i) - xwall_phys(i-1);
    p_avg = (p_wall(i) + p_wall(i-1)) / 2;
    F_vertical = F_vertical + p_avg * dx * w;
end

% Also need settling chamber force contribution
% In settling chamber: p_sc ≈ p0 (since M_sc is very small)
% The settling chamber contributes to axial thrust as well
p_sc = p0_force * (1 + (G-1)/2 * M_sc^2)^(-G/(G-1));  % ≈ p0

fprintf('  Vertical force on one nozzle wall (expansion section) = %.2f N\n', F_vertical);
fprintf('  (This is the force due to pressure on a single nozzle plane)\n');

%% FUNCTIONS
function [M, Nu, Mu] = PMF(G, M_in, Nu_in, ~)
    if M_in >= 1
        nu_rad = sqrt((G+1)/(G-1)) * atan(sqrt((G-1)/(G+1)*(M_in^2-1))) - atan(sqrt(M_in^2-1));
        Nu = rad2deg(nu_rad);
        Mu = rad2deg(asin(1/M_in));
        M = M_in;
    elseif Nu_in > 0
        % Inverse: solve for M given Nu using fzero
        fun = @(Mach) rad2deg(sqrt((G+1)/(G-1)) * atan(sqrt((G-1)/(G+1)*(Mach^2-1))) - atan(sqrt(Mach^2-1))) - Nu_in;
        M = fzero(fun, [1.0001 30]);
        Nu = Nu_in;
        Mu = rad2deg(asin(1/M));
    else
        M = 1; Nu = 0; Mu = 90;
    end
end

function [xarc, yarc] = Arc(R, ThetaArc)
    theta = deg2rad(ThetaArc);
    xarc = R * sin(theta);
    yarc = R * (1 - cos(theta));
end