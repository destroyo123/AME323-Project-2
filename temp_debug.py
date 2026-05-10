import math

gamma = 1.4
M_e = 6.0
nu = lambda M: math.sqrt((gamma+1)/(gamma-1))*math.degrees(math.atan(math.sqrt(((gamma-1)/(gamma+1))*(M*M-1)))) - math.degrees(math.atan(math.sqrt(M*M-1)))
nu_e = nu(M_e)
delta_max = nu_e/2
N = 40
ddelta = delta_max/N
h_t = 0.0033
xThroat = [h_t*math.sin(i*ddelta*math.pi/180) for i in range(N+1)]
yThroat = [h_t*(2-math.cos(i*ddelta*math.pi/180)) for i in range(N+1)]
nu0 = [(i+1)*ddelta for i in range(N)]
theta0 = list(nu0)
Kp0 = [nu0[i]+theta0[i] for i in range(N)]
Km0 = [nu0[i]-theta0[i] for i in range(N)]

# Use a numeric invert to solve M from nu

def M_from_nu(nu_val):
    lo, hi = 1.0, 1000.0
    for _ in range(120):
        mid = 0.5*(lo+hi)
        val = math.sqrt((gamma+1)/(gamma-1))*math.degrees(math.atan(math.sqrt(((gamma-1)/(gamma+1))*(mid*mid-1)))) - math.degrees(math.atan(math.sqrt(mid*mid-1)))
        if val > nu_val:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo+hi)

mu0 = [math.degrees(math.asin(1.0/M_from_nu(n))) for n in nu0]

x0 = xThroat[1:]
y0 = yThroat[1:]

x_center = []
for i in range(N):
    slope_down = math.tan(math.radians(theta0[i]-mu0[i]))
    if slope_down >= 0:
        slope_down = math.tan(math.radians(theta0[i]+mu0[i]))
        slope_down = -abs(slope_down)
    x_c = x0[i] + (-y0[i]/slope_down)
    x_center.append(x_c)

Kp_center = list(Kp0)
Km_center = list(Kp_center)

x_trans = [[None]*N for _ in range(N)]
y_trans = [[None]*N for _ in range(N)]
M_trans = [[None]*N for _ in range(N)]
nu_trans = [[None]*N for _ in range(N)]
delta_trans = [[None]*N for _ in range(N)]
mu_trans = [[None]*N for _ in range(N)]

for j in range(1,N):
    for i in range(j):
        Kp = Kp0[j]
        Km = Km0[i]
        nu_t = 0.5*(Kp+Km)
        delta_t = 0.5*(Kp-Km)
        M_t = M_from_nu(nu_t)
        mu_t = math.degrees(math.asin(1.0/M_t))
        nu_trans[i][j] = nu_t
        delta_trans[i][j] = delta_t
        M_trans[i][j] = M_t
        mu_trans[i][j] = mu_t
        if i == 0:
            xA, yA = x_center[j], 0.0
            delta_A, mu_A = 0.0, math.degrees(math.asin(1.0/M_from_nu(Kp_center[j])))
        else:
            xA, yA = x_trans[i-1][j], y_trans[i-1][j]
            delta_A = delta_trans[i-1][j]
            mu_A = mu_trans[i-1][j]
        if j == i+1:
            xB, yB = x_center[i], 0.0
            delta_B, mu_B = 0.0, math.degrees(math.asin(1.0/M_from_nu(Kp_center[i])))
        else:
            xB, yB = x_trans[i][j-1], y_trans[i][j-1]
            delta_B = delta_trans[i][j-1]
            mu_B = mu_trans[i][j-1]
        slope_Cplus = math.tan(math.radians(delta_A - mu_A))
        slope_Cminus = math.tan(math.radians(delta_B + mu_B))
        if abs(slope_Cplus - slope_Cminus) < 1e-12:
            x_int = 0.5*(xA + xB)
        else:
            x_int = (yB - yA + slope_Cplus*xA - slope_Cminus*xB)/(slope_Cplus - slope_Cminus)
        y_int = yA + slope_Cplus*(x_int - xA)
        x_trans[i][j] = x_int
        y_trans[i][j] = y_int

print('x_center last', x_center[-1])
print('x_trans sample')
for i in [0,1,2,3,10,20]:
    for j in [i+1, min(i+2,N-1), min(i+5,N-1)]:
        print('i', i+1,'j',j+1,'x', x_trans[i][j])

x_wall_start = xThroat[-1]
y_wall_start = yThroat[-1]
delta_wall_start = delta_max
x_wall_pts = [x_wall_start]
y_wall_pts = [y_wall_start]
delta_wall_pts = [delta_wall_start]
M_wall=[]
nu_wall=[]
mu_wall=[]

for k in range(1,N+1):
    i = N - k
    Km_w = Km_center[i]
    Kp_w = Kp0[-1]
    nu_w = 0.5*(Kp_w + Km_w)
    delta_w = 0.5*(Kp_w - Km_w)
    nu_wall.append(nu_w)
    delta_wall.append(delta_w)
    M_w = M_from_nu(nu_w)
    mu_w = math.degrees(math.asin(1.0/M_w))
    M_wall.append(M_w)
    mu_wall.append(mu_w)
    if i < N-1 and x_trans[i][N-1] is not None:
        xB, yB = x_trans[i][N-1], y_trans[i][N-1]
        delta_B, mu_B = delta_trans[i][N-1], mu_trans[i][N-1]
    else:
        xB, yB = x_center[i], 0.0
        delta_B, mu_B = 0.0, math.degrees(math.asin(1.0/M_from_nu(Kp_center[i])))
    slope_Cm = 0.5*(math.tan(math.radians(delta_B + mu_B)) + math.tan(math.radians(delta_w + mu_w)))
    slope_wall = math.tan(math.radians(0.5*(delta_wall_pts[k-1] + delta_w)))
    if abs(slope_Cm - slope_wall) < 1e-12:
        x_w = x_wall_pts[k-1] + 0.01*h_t
    else:
        x_w = (yB - y_wall_pts[k-1] + slope_wall*x_wall_pts[k-1] - slope_Cm*xB)/(slope_wall - slope_Cm)
    y_w = yB + slope_Cm*(x_w - xB)
    if x_w <= x_wall_pts[k-1]:
        x_w = x_wall_pts[k-1] + 0.01*h_t
        y_w = y_wall_pts[k-1] + slope_wall*0.01*h_t
    x_wall_pts.append(x_w)
    y_wall_pts.append(y_w)
    delta_wall_pts.append(delta_w)
    if x_w > 1000:
        print('large wall',k,i+1,x_w,y_w,delta_w)

print('wall points', len(x_wall_pts))
print('first few wall', list(zip(x_wall_pts[:5], y_wall_pts[:5], delta_wall_pts[:5])))
print('last few wall', list(zip(x_wall_pts[-5:], y_wall_pts[-5:], delta_wall_pts[-5:])))
