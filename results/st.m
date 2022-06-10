% Stability calculations
dt = 1;
zeta = 1.9; % zeta > 1


% mine w_i, alpha =1
w_i = 1/dt * ( log(2*zeta) / log(exp(-0.5) + exp(0.5)) )

% taflove w_i, alpha=1
w_t = -2/dt * log(zeta + sqrt(zeta^2 + 1))

MU_0 = 1.2566371e-6;
EPS_0 = 8.85418781762e-12;
C_CONST = 299792458.0;
dz = 0.02e-6;
alpha = 0.98;
dt = 2^(1-1/alpha) * (sqrt(EPS_0*MU_0) * dz)^(1/alpha)
classic_dt = sqrt(EPS_0*MU_0)*dz