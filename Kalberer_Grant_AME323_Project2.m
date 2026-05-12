%% AME 323 Project 2 - Mach 6 Blowdown Wind Tunnel / 2-D Nozzle Design
% Self-contained MATLAB script. Run this file directly.
%
% This script completes the computational prompts in Project 2:
%   1) method-of-characteristics (MOC) top-half nozzle contour design,
%      with a smooth circular-arc throat of radius equal to throat half-height;
%   2) grid convergence, nozzle/flow-field plots, centerline distributions,
%      required pressure ratio, and uniform-flow triangle length;
%   3) tunnel operating pressure, mass flow, run time, unit Reynolds number,
%      settling chamber size, and nozzle force estimates.
%
% Coordinate system: x starts at the sonic throat. y=0 is the centerline.
% Only the top half of the nozzle is generated and plotted.

clear; close all; clc;

%% ----------------------- User-adjustable inputs ------------------------
gamma  = 1.40;                 % air specific heat ratio [-]
Rgas   = 287.05287;            % air gas constant [J/(kg K)]
Me     = 6.0;                  % design exit Mach number [-]
N      = 65;                   % final design number of initial PM waves [-]
Nstudy = [8 12 16 24 32 48 65 80 100 130]; % convergence study wave counts

exitHeight = 0.350;            % full nozzle exit height [m]
width      = 0.500;            % nozzle/tunnel width [m]
pb_psia    = 0.100;            % maintained test-section static pressure [psia]
p0max_psia = 500.0;            % maximum regulator total pressure [psia]
Tmin_free  = 55.0;             % minimum allowable free-stream static T [K]

tankVolume = 40.0;             % compressed-air tank volume [m^3]
tankTemp   = 300.0;            % approximately constant tank temperature [K]
tankInitialPressure_psia = p0max_psia; % assumed full tank pressure [psia]

figDir = fullfile(pwd,'Project2_figures');
if ~exist(figDir,'dir'), mkdir(figDir); end

%% ---------------------- Basic design quantities ------------------------
nu_e       = prandtlMeyer(Me,gamma);
theta_max  = 0.5*nu_e;                 % maximum wall angle [rad]
AR_e       = areaMach(Me,gamma);        % A_e/A* = h_e/h_t for constant width
hExit      = 0.5*exitHeight;            % top-half exit height [m]
hThroat    = hExit/AR_e;                % top-half throat height [m]
Rthroat    = hThroat;                   % required circular throat radius [m]
Ae         = exitHeight*width;          % full exit area [m^2]
Athroat    = Ae/AR_e;                   % full throat area [m^2]

T0min      = Tmin_free*(1+(gamma-1)/2*Me^2);
pp0_exit   = isenPressureRatio(Me,gamma);      % p_e/p_0
TT0_exit   = isenTemperatureRatio(Me,gamma);   % T_e/T_0
p0_over_pe = 1/pp0_exit;
p0min_psia = pb_psia/pp0_exit;

fprintf('\nAME 323 Project 2 - Mach %.1f nozzle design\n',Me);
fprintf('Exit area ratio A_e/A*                         = %.6f\n',AR_e);
fprintf('Exit full height x width                         = %.3f m x %.3f m\n',exitHeight,width);
fprintf('Throat full height                               = %.6f m\n',2*hThroat);
fprintf('Throat circular-arc radius                       = %.6f m\n',Rthroat);
fprintf('Maximum wall angle                               = %.3f deg\n',rad2deg(theta_max));
fprintf('Minimum total temperature for T_e >= %.1f K       = %.3f K\n',Tmin_free,T0min);
fprintf('p_e/p_0 at Mach %.1f                              = %.8f\n',Me,pp0_exit);
fprintf('p_0/p_e required for ideal expansion              = %.3f\n',p0_over_pe);
fprintf('Minimum p_0 for p_b = %.3f psia                  = %.3f psia\n',pb_psia,p0min_psia);

%% ------------------------- Grid convergence ----------------------------
err_AR = zeros(size(Nstudy));
Lstudy = zeros(size(Nstudy));
for k = 1:numel(Nstudy)
    Dk = designNozzleMOC(Nstudy(k),Me,gamma,hThroat);
    err_AR(k) = abs(Dk.areaRatioGeom - AR_e)/AR_e;
    Lstudy(k) = Dk.xExit;
end

figure('Name','Grid convergence');
loglog(Nstudy,err_AR,'o-','LineWidth',1.4); grid on;
xlabel('Number of initial Prandtl-Meyer waves, N');
ylabel('| (h_e/h_t)_{MOC} - (A_e/A*)_{isentropic} | / (A_e/A*)');
title('Grid convergence of geometric exit area ratio');
saveas(gcf,fullfile(figDir,'01_grid_convergence.png'));

%% ------------------------- Final MOC design ----------------------------
D = designNozzleMOC(N,Me,gamma,hThroat);

fprintf('\nFinal design using N = %d initial waves\n',N);
fprintf('MOC geometric exit area ratio h_e/h_t             = %.6f\n',D.areaRatioGeom);
fprintf('Relative area-ratio error                         = %.4e\n',abs(D.areaRatioGeom-AR_e)/AR_e);
fprintf('Nozzle expansion length                           = %.6f m\n',D.xExit);
fprintf('Uniform-flow triangle length inside nozzle        = %.6f m\n',D.uniformTriangleLength);

%% --------------------------- Nozzle plots ------------------------------
% Top-half smooth wall contour only; circular throat arc plus MOC wall.
figure('Name','Top-half nozzle contour'); hold on; grid on; axis equal;
plot(D.wall.x,D.wall.y,'k-','LineWidth',2.0);
plot([0 D.xExit],[0 0],'k--','LineWidth',1.0);
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Smooth top-half Mach %.1f nozzle contour, N = %d',Me,N));
legend('Top wall contour','Centerline','Location','best');
saveas(gcf,fullfile(figDir,'02_top_half_nozzle_contour.png'));

% Mach-number contour map in top half.
plotScalarField(D,'M','Mach number',fullfile(figDir,'03_mach_contours_top_half.png'));
plotScalarField(D,'p_p0','p/p_0',fullfile(figDir,'04_pressure_ratio_contours_top_half.png'));
plotScalarField(D,'T_T0','T/T_0',fullfile(figDir,'05_temperature_ratio_contours_top_half.png'));

% Centerline distributions.
figure('Name','Centerline distributions');
plot(D.centerline.x,D.centerline.M,'o-','LineWidth',1.4); hold on; grid on;
plot(D.centerline.x,D.centerline.p_p0,'s-','LineWidth',1.4);
plot(D.centerline.x,D.centerline.T_T0,'^-','LineWidth',1.4);
xlabel('x [m]'); ylabel('Nondimensional value');
title('Centerline distributions from throat to exit');
legend('M','p/p_0','T/T_0','Location','best');
saveas(gcf,fullfile(figDir,'06_centerline_distributions.png'));

%% --------------------- Tunnel operating calculations -------------------
psia2Pa = 6894.757293168;
p0vec_psia = linspace(max(p0min_psia,1e-6),p0max_psia,250);
p0vec = p0vec_psia*psia2Pa;
T0 = T0min;

mdot = massFlowChoked(p0vec,T0,Athroat,gamma,Rgas);

pTank0 = tankInitialPressure_psia*psia2Pa;
% Available tank mass above the selected regulated pressure. This assumes the
% run ends when tank pressure falls to the regulator set pressure.
mAvailable = max((pTank0 - p0vec)*tankVolume/(Rgas*tankTemp),0);
runTime = mAvailable ./ mdot;

Te = T0*TT0_exit;
pe = p0vec*pp0_exit;
rhoe = pe./(Rgas*Te);
Ue = Me*sqrt(gamma*Rgas*Te);
mue = sutherlandMu(Te);
unitRe = rhoe.*Ue./mue;

Msc = 0.05;
Asc = Athroat*areaMach(Msc,gamma);
hsc = Asc/width;
Tsc = T0*isenTemperatureRatio(Msc,gamma);
psc_over_p0 = isenPressureRatio(Msc,gamma);
Usc = Msc*sqrt(gamma*Rgas*Tsc);

fprintf('\nTunnel characteristics\n');
fprintf('Minimum total pressure for shock-free ideal exit  = %.3f psia\n',p0min_psia);
fprintf('Operating total-pressure range used in plots      = %.3f to %.3f psia\n',max(p0min_psia,1e-6),p0max_psia);
fprintf('Mass flow at p0_min and 500 psia                  = %.3f, %.3f kg/s\n',mdot(1),mdot(end));
fprintf('Run time at p0_min and 500 psia                   = %.3f, %.3f s\n',runTime(1),runTime(end));
fprintf('Free-stream unit Re at p0_min and 500 psia        = %.3e, %.3e 1/m\n',unitRe(1),unitRe(end));
fprintf('Settling chamber area for M <= 0.05               = %.6f m^2\n',Asc);
fprintf('Settling chamber height for width %.3f m          = %.6f m\n',width,hsc);

figure('Name','Mass flow range');
plot(p0vec_psia,mdot,'LineWidth',1.5); grid on;
xlabel('Regulated total pressure p_0 [psia]'); ylabel('\dot{m} [kg/s]');
title(sprintf('Choked mass flow, T_0 = %.1f K',T0));
saveas(gcf,fullfile(figDir,'07_mass_flow_vs_p0.png'));

figure('Name','Run time');
plot(p0vec_psia,runTime,'LineWidth',1.5); grid on;
xlabel('Regulated total pressure p_0 [psia]'); ylabel('Estimated run time [s]');
title(sprintf('Run time with %.1f m^3 tank initially at %.0f psia, %.0f K',tankVolume,tankInitialPressure_psia,tankTemp));
saveas(gcf,fullfile(figDir,'08_runtime_vs_p0.png'));

figure('Name','Unit Reynolds number');
plot(p0vec_psia,unitRe,'LineWidth',1.5); grid on;
xlabel('Regulated total pressure p_0 [psia]'); ylabel('Re'' [1/m]');
title(sprintf('Free-stream unit Reynolds number, M = %.1f, T_e = %.1f K',Me,Te));
saveas(gcf,fullfile(figDir,'09_unit_re_vs_p0.png'));

%% ----------------------- Nozzle force estimates ------------------------
% Momentum/thrust contribution from the settling chamber to the exit.
Fx_momentum = mdot.*(Ue - Usc);

% Pressure force normal to one top nozzle wall in the designed expansion.
% This uses gauge pressure relative to the maintained test-section pressure.
[Fy_top, Fx_pressure_wall] = wallPressureForces(D,p0vec,pb_psia*psia2Pa,width,gamma);

figure('Name','Nozzle forces');
plot(p0vec_psia,Fx_momentum,'LineWidth',1.5); hold on; grid on;
plot(p0vec_psia,Fy_top,'LineWidth',1.5);
plot(p0vec_psia,Fx_pressure_wall,'LineWidth',1.5);
xlabel('Regulated total pressure p_0 [psia]'); ylabel('Force [N]');
title('Estimated nozzle force components');
legend('Axial momentum force, settling chamber to exit',...
       'Vertical pressure force on one top expansion wall',...
       'Axial pressure force on one top expansion wall','Location','best');
saveas(gcf,fullfile(figDir,'10_force_estimates.png'));

fprintf('\nNozzle force estimates at p0 = %.1f psia\n',p0vec_psia(end));
fprintf('Axial momentum force, settling chamber to exit    = %.3f N\n',Fx_momentum(end));
fprintf('Vertical pressure force on one top expansion wall = %.3f N (downward on top wall)\n',Fy_top(end));
fprintf('Axial pressure force on one top expansion wall    = %.3f N\n',Fx_pressure_wall(end));

%% ----------------------------- Summary ---------------------------------
fprintf('\nFigures written to: %s\n',figDir);
fprintf('Done. Use the printed values and PNG figures in the design report.\n');

%% =======================================================================
%                               FUNCTIONS
% ========================================================================

function D = designNozzleMOC(N,Me,g,hT)
%DESIGNNOZZLEMOC Top-half 2-D MOC nozzle with a circular-arc throat.
% The initial expansion is represented by N equally spaced PM waves on a
% circular throat arc of radius hT. The straightening wall is generated by
% compatibility with a uniform Mach-Me exit stream.

    nu_e = prandtlMeyer(Me,g);
    thMax = 0.5*nu_e;

    % Avoid exactly sonic throat, where mu = 90 deg. Include a tiny first
    % wave so the wall starts smoothly just downstream of the throat.
    th = linspace(thMax/N,thMax,N);   % wall angle = PM angle on simple wave
    nu = th;
    M  = invPrandtlMeyer(nu,g);
    mu = asin(1./M);

    % Smooth circular throat arc. At x=0,y=hT the wall is horizontal.
    W.x = hT*sin(th);
    W.y = 2*hT - hT*cos(th);  % = hT*(2 - cos(th))
    W.theta = th;
    W.nu = nu;
    W.M = M;
    W.mu = mu;

    % Triangular MOC grid. P(i,j) is intersection of C+ from centerline
    % wave i with C- wave j, for 1 <= i <= j <= N. Diagonal is centerline.
    P = repmat(struct('x',NaN,'y',NaN,'theta',NaN,'nu',NaN,'M',NaN,'mu',NaN,'Kp',NaN,'Km',NaN),N,N);

    for j = 1:N
        for i = 1:j
            Kp = -2*th(i);     % C+ invariant from centerline/reflection
            Km =  2*th(j);     % C- invariant from circular-arc expansion
            theta_ij = 0.5*(Kp + Km);
            nu_ij    = 0.5*(Km - Kp);
            M_ij     = invPrandtlMeyer(nu_ij,g);
            mu_ij    = asin(1/M_ij);

            if i == 1
                srcMinus = WPoint(W,j);
            else
                srcMinus = P(i-1,j);
            end

            if i == j
                % Intersect the incoming C- segment with the centerline y=0.
                mMinus = tan(0.5*((srcMinus.theta - srcMinus.mu) + (theta_ij - mu_ij)));
                x_ij = srcMinus.x - srcMinus.y/mMinus;
                y_ij = 0;
            else
                srcPlus = P(i,j-1);
                mMinus = tan(0.5*((srcMinus.theta - srcMinus.mu) + (theta_ij - mu_ij)));
                mPlus  = tan(0.5*((srcPlus.theta  + srcPlus.mu)  + (theta_ij + mu_ij)));
                [x_ij,y_ij] = intersectLines(srcMinus.x,srcMinus.y,mMinus,srcPlus.x,srcPlus.y,mPlus);
            end

            P(i,j).x = x_ij; P(i,j).y = y_ij;
            P(i,j).theta = theta_ij; P(i,j).nu = nu_ij;
            P(i,j).M = M_ij; P(i,j).mu = mu_ij;
            P(i,j).Kp = Kp; P(i,j).Km = Km;
        end
    end

    % Straightening wall. Each wall point is reached by a C+ from P(i,N)
    % and is compatible with the uniform exit stream (theta=0, nu=nu_e).
    Q = repmat(struct('x',NaN,'y',NaN,'theta',NaN,'nu',NaN,'M',NaN,'mu',NaN),1,N);
    prev.x = W.x(end); prev.y = W.y(end); prev.theta = W.theta(end);
    for i = 1:N
        src = P(i,N);
        Kp = src.Kp;
        theta_q = 0.5*(nu_e + Kp);
        nu_q    = 0.5*(nu_e - Kp);
        M_q     = invPrandtlMeyer(nu_q,g);
        mu_q    = asin(1/M_q);

        mChar = tan(0.5*((src.theta + src.mu) + (theta_q + mu_q)));
        mWall = tan(0.5*(prev.theta + theta_q));
        [x_q,y_q] = intersectLines(src.x,src.y,mChar,prev.x,prev.y,mWall);

        Q(i).x = x_q; Q(i).y = y_q; Q(i).theta = theta_q;
        Q(i).nu = nu_q; Q(i).M = M_q; Q(i).mu = mu_q;
        prev = Q(i);
    end

    % Combine wall coordinates; include exact throat point.
    wallX = [0 W.x Q.x];
    wallY = [hT W.y Q.y];
    wallTheta = [0 W.theta Q.theta];
    wallM = [1 W.M Q.M];

    % Grid point clouds for contouring. Include centerline, interior, wall.
    x = []; y = []; Mvals = []; thetas = [];
    x(end+1) = 0; y(end+1) = 0; Mvals(end+1) = 1; thetas(end+1) = 0;
    for j = 1:N
        for i = 1:j
            x(end+1) = P(i,j).x; %#ok<AGROW>
            y(end+1) = P(i,j).y; %#ok<AGROW>
            Mvals(end+1) = P(i,j).M; %#ok<AGROW>
            thetas(end+1) = P(i,j).theta; %#ok<AGROW>
        end
    end
    x = [x wallX]; y = [y wallY]; Mvals = [Mvals wallM]; thetas = [thetas wallTheta];

    p_p0 = isenPressureRatio(Mvals,g);
    T_T0 = isenTemperatureRatio(Mvals,g);

    % Exit coordinates.
    xExit = Q(end).x; yExit = Q(end).y;

    % Centerline values from throat to exit. After the final centerline
    % characteristic, the flow is uniform, so M, p/p0, and T/T0 remain fixed
    % until the geometric exit plane.
    cx = zeros(1,N+2); cM = zeros(1,N+2);
    cx(1) = 0; cM(1) = 1;
    for j = 1:N
        cx(j+1) = P(j,j).x;
        cM(j+1) = P(j,j).M;
    end
    cx(end) = xExit; cM(end) = Me;

    % Length of triangular uniform-flow region: intersection of final C+
    % from centerline with the straight exit wall, measured from the exit.
    finalCenter = P(N,N);
    mFinalPlus = tan(finalCenter.theta + finalCenter.mu);
    xHitExitHeight = finalCenter.x + (yExit - finalCenter.y)/mFinalPlus;
    uniformTriangleLength = max(xExit - xHitExitHeight,0);

    D.P = P; D.W = W; D.Q = Q;
    D.wall.x = wallX; D.wall.y = wallY; D.wall.theta = wallTheta; D.wall.M = wallM;
    D.cloud.x = x; D.cloud.y = y; D.cloud.M = Mvals; D.cloud.theta = thetas;
    D.cloud.p_p0 = p_p0; D.cloud.T_T0 = T_T0;
    D.centerline.x = cx; D.centerline.M = cM;
    D.centerline.p_p0 = isenPressureRatio(cM,g);
    D.centerline.T_T0 = isenTemperatureRatio(cM,g);
    D.xExit = xExit; D.yExit = yExit; D.hThroat = hT;
    D.areaRatioGeom = yExit/hT;
    D.uniformTriangleLength = uniformTriangleLength;
end

function S = WPoint(W,j)
    S.x = W.x(j); S.y = W.y(j); S.theta = W.theta(j); S.nu = W.nu(j); S.M = W.M(j); S.mu = W.mu(j);
end

function [x,y] = intersectLines(x1,y1,m1,x2,y2,m2)
    % y-y1=m1(x-x1), y-y2=m2(x-x2)
    den = (m1 - m2);
    if abs(den) < 1e-13
        x = 0.5*(x1+x2);
    else
        x = (y2 - y1 + m1*x1 - m2*x2)/den;
    end
    y = y1 + m1*(x-x1);
end

function nu = prandtlMeyer(M,g)
    nu = sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(M.^2-1))) - atan(sqrt(M.^2-1));
end

function M = invPrandtlMeyer(nu,g)
    % Robust bisection inversion of the PM function. nu can be scalar/vector.
    M = zeros(size(nu));
    for k = 1:numel(nu)
        if nu(k) <= 0
            M(k) = 1;
        else
            lo = 1 + 1e-10; hi = 50;
            while prandtlMeyer(hi,g) < nu(k)
                hi = 2*hi;
            end
            for it = 1:80
                mid = 0.5*(lo+hi);
                if prandtlMeyer(mid,g) < nu(k)
                    lo = mid;
                else
                    hi = mid;
                end
            end
            M(k) = 0.5*(lo+hi);
        end
    end
end

function A = areaMach(M,g)
    A = (1./M).*((2/(g+1))*(1+(g-1)/2*M.^2)).^((g+1)/(2*(g-1)));
end

function pp0 = isenPressureRatio(M,g)
    pp0 = (1+(g-1)/2*M.^2).^(-g/(g-1));
end

function TT0 = isenTemperatureRatio(M,g)
    TT0 = 1./(1+(g-1)/2*M.^2);
end

function mdot = massFlowChoked(p0,T0,Astar,g,R)
    mdot = Astar.*p0./sqrt(T0).*sqrt(g/R).*(2/(g+1))^((g+1)/(2*(g-1)));
end

function mu = sutherlandMu(T)
    muRef = 1.716e-5; Tref = 273.15; S = 110.4;
    mu = muRef*(T./Tref).^(3/2).*((Tref+S)./(T+S));
end

function plotScalarField(D,fieldName,labelText,fileName)
    x = D.cloud.x(:); y = D.cloud.y(:); z = D.cloud.(fieldName)(:);

    % Remove any duplicate/invalid points before interpolation.
    ok = isfinite(x) & isfinite(y) & isfinite(z);
    x = x(ok); y = y(ok); z = z(ok);
    [xyUnique,ia] = unique(round([x y],10),'rows','stable'); %#ok<ASGLU>
    x = x(ia); y = y(ia); z = z(ia);

    xq = linspace(min(x),max(x),280);
    yq = linspace(0,max(y),140);
    [Xq,Yq] = meshgrid(xq,yq);
    Zq = griddata(x,y,z,Xq,Yq,'natural');

    % Mask points above the smooth top wall.
    yWall = interp1(D.wall.x,D.wall.y,xq,'linear','extrap');
    for k = 1:numel(xq)
        Zq(Yq(:,k) > yWall(k),k) = NaN;
    end

    figure('Name',labelText); hold on;
    contourf(Xq,Yq,Zq,30,'LineColor','none'); colorbar;
    plot(D.wall.x,D.wall.y,'k-','LineWidth',2.0);
    plot([0 D.xExit],[0 0],'k--','LineWidth',1.0);
    axis equal tight; grid on;
    xlabel('x [m]'); ylabel('y [m]'); title([labelText ' in top half of nozzle']);
    saveas(gcf,fileName);
end

function [Fy_top,Fx_wall] = wallPressureForces(D,p0vec,pb,width,g)
    x = D.wall.x(:); y = D.wall.y(:); M = D.wall.M(:);
    % Use gauge pressure relative to external/test-section pressure.
    pp0 = isenPressureRatio(M,g);
    dx = diff(x); dy = diff(y);
    Mmid_pp0 = 0.5*(pp0(1:end-1)+pp0(2:end));
    slope = dy./max(dx,eps);

    Fy_top = zeros(size(p0vec));
    Fx_wall = zeros(size(p0vec));
    for k = 1:numel(p0vec)
        pg = p0vec(k)*Mmid_pp0 - pb;
        % For a top wall, the pressure force from the flow on the wall is
        % downward and generally downstream/opposing the outward normal. The
        % magnitudes below are useful for structural load estimates.
        Fy_top(k) = width*sum(pg.*dx);          % vertical magnitude [N]
        Fx_wall(k) = width*sum(pg.*dy);         % axial component magnitude [N]
    end
end
