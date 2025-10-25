function [XDOT] = Nonlinear_6DOF_Model(X, U)

% State vector
x1 = X(1); % u (body-x velocity)
x2 = X(2); % v (body-y velocity)
x3 = X(3); % w (body-z velocity)
x4 = X(4); % p (roll rate)
x5 = X(5); % q (pitch rate)
x6 = X(6); % r (yaw rate)
x7 = X(7); % phi (roll angle)
x8 = X(8); % theta (pitch angle)
x9 = X(9); % psi (yaw angle)

u1 = U(1); % Aileron deflection (delta_a)
u2 = U(2); % Elevator deflection (delta_e)
u3 = U(3); % Rudder deflection (delta_r)
u4 = U(4); % Engine 1 throttle (tau1)
u5 = U(5); % Engine 2 throttle (tau2)

%% Constants
m = 120000;
cbar = 6.6;
lt = 24.8;
S = 260;
St = 64;

Xcg = 0.23 * cbar;
Ycg = 0;
Zcg = 0.10 * cbar;

Xac = 0.12 * cbar;
Yac = 0;
Zac = 0;

Xapt1 = 0;
Yapt1 = -7.94;
Zapt1 = -1.9;

Xapt2 = 0;
Yapt2 = 7.94;
Zapt2 = -1.9;

%% Other constants
rho = 1.225;
g = 9.81;
depsda = 0.25;

a3 = -768.5;
a2 = 609.2;
a1 = -155.2;
a0 = 15.212;

alpha_LO = -11.5 * (pi/180);
alpha_switch = 14.5 * (pi/180);
n = 5.5;

%% Intermediate Variables
% Calculated airspeed
va = sqrt(x1^2 + x2^2 + x3^2);

va_min = 1e-3; % 0.001 m/s
va_for_division = max(va, va_min);
% ------------------------------------------------------------------

% Calculate alpha and beta
alpha = atan2(x3, x1);
% Use clamped velocity for beta calculation
beta = asin(x2 / va_for_division);

% Calculate dynamic pressure
Q = 0.5 * rho * va^2;

%  define the vector omega and Vb
omega_b = [x4; x5; x6];
Vb = [x1; x2; x3];

%% Aerodynamic Force Coefficients
if alpha <= alpha_switch
    CL_wb = n * (alpha - alpha_LO);
else
    CL_wb = a3 * alpha^3 + a2 * alpha^2 + a1 * alpha + a0;
end

% Calculate lift coefficient for tail
epsilon = depsda * (alpha - alpha_LO);
% Use clamped velocity for the dynamic term (lt / va) in alpha_t
alpha_t = alpha - epsilon + u2 + 1.3 * x5 * (lt / va_for_division); 
CL_t = 3.1 * (St / S) * alpha_t;

% Total lift force
CL = CL_wb + CL_t;

% Total Drag Force
CD = 0.13 + 0.07 * (5.5 * alpha + 0.654)^2;

% Calculate side force
CY = -1.6 * beta + 0.24 * u3;

%% Dimensional Aerodynamic Forces
FA_s = [-CD * Q * S;
         CY * Q * S;
        -CL * Q * S];

% Rotate the force to body axis
C_bs = [cos(alpha)  0  -sin(alpha);
         0          1   0;
         sin(alpha) 0   cos(alpha)];

FA_b = C_bs * FA_s;

%% Aerodynamic moments about the AC
eta11 = -1.4 * beta;
eta21 = -0.59 - (3.1 * (St * lt) / (S * cbar)) * (alpha - epsilon);
eta31 = (1 - alpha * (180 / (15 * pi))) * beta;

eta = [eta11;
        eta21;
        eta31];

% dCMdx matrix uses va in the original formulation: (cbar / va)
dCMdx = (cbar / va_for_division) * [ -11   0    5;
                         0   (-4.03 * (St * lt^2) / (S * cbar^2))   0;
                         1.7   0   -11.5 ];

dCMdu = [ -0.6   0     0.22;
           0   (-3.1 * (St * lt) / (S * cbar))   0;
           0     0    -0.63 ];

% Calculate CM about aerodynamic center in Fb
CMac_b = eta + dCMdx * omega_b + dCMdu * [u1; u2; u3];

%% Aerodynamic Moments about AC
MAac_b = CMac_b * Q * S * cbar;

%% Aerodynamic Moments about CG
rcg_b = [Xcg; Ycg; Zcg];
rac_b = [Xac; Yac; Zac];

MAcg_b = MAac_b + cross(FA_b, rcg_b - rac_b);

%% Engine Forces and Moments
% u4 and u5 are throttle settings (tau), not forces. F = tau * m * g
F1 = u4 * m * g; 
F2 = u5 * m * g;

FE1_b = [F1; 0; 0];
FE2_b = [F2; 0; 0];
FE_b = FE1_b + FE2_b;

% Engine moment due to offset of engine thrust from CoG
mew1 = [Xcg - Xapt1;
        Yapt1 - Ycg;
        Zcg - Zapt1];

mew2 = [Xcg - Xapt2;
        Yapt2 - Ycg;
        Zcg - Zapt2];

MEcg1_b = cross(mew1, FE1_b);
MEcg2_b = cross(mew2, FE2_b);
MEcg_b = MEcg1_b + MEcg2_b;

%% Gravity Effects
g_b = [-g*sin(x8);
       g*cos(x8)*sin(x7);
       g*cos(x8)*cos(x7)];

Fg_b = m * g_b;

%% Inertia Matrix
Ib_m = [40.07      0       -2.0923;
         0        64        0;
        -2.0923   0       99.92];

invIb = (1/m) * [0.0249836      0       0.000523151;
                   0       0.015625      0;
              0.000523151     0      0.010019];

%% Translational Dynamics
F_b = Fg_b + FE_b + FA_b;
x1to3dot = (1/m) * F_b - cross(omega_b, Vb);

%% Rotational Dynamics
Mcg_b = MAcg_b + MEcg_b;
x4to6dot = invIb * (Mcg_b - cross(omega_b, Ib_m * omega_b));

%% Kinematics (Euler angles)
H_phi = [1              sin(x7) * tan(x8)      cos(x7) * tan(x8);
         0              cos(x7)               -sin(x7);
         0              sin(x7)/cos(x8)       cos(x7)/cos(x8)];

x7to9dot = H_phi * omega_b;

%% Assemble state derivative
XDOT = [x1to3dot;
        x4to6dot;
        x7to9dot];

end