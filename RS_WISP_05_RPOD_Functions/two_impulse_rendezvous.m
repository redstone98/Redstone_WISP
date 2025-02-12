function [Delta_v_0, Delta_v_f, delta_r_t_mat, delta_v_t_mat, v_plus_ECI, t_vector] = two_impulse_rendezvous(r_target, v_target, r_chaser, v_chaser, tf)

r_0 = r_target;
v_0 = v_target;
r = r_chaser;
v = v_chaser;

% Reference Frame of the Target
i_hat = r_0/norm(r_0);
j_hat = v_0/norm(v_0);
k_hat = cross(i_hat, j_hat);

% Transformation Matrix from ECI to Space Station Frame
Q_LVLH_ECI = [i_hat';j_hat';k_hat'];

% Position vector of the spacecraft relative to the space station (ECI)
delta_r = r - r_0;
n = norm(v_0)/norm(r_0);
Omega_target = n * k_hat;
delta_v = v - v_0 - cross(Omega_target, delta_r);

% Relative position vector at the beginning of the rendezvous maneuver
delta_r_0 = Q_LVLH_ECI * delta_r;

% Relative velocity just before launch into the rendezvous trajectory is
delta_v_0_minus = Q_LVLH_ECI * delta_v;

% Calculate Clohessy-Whiltshire matrix for t = t_f = 28800s and n
nt = n*tf;

Phi_rr_tf = [4 - 3*cos(nt),   0, 0;
             6*(sin(nt) - nt),1,0;
             0,0, cos(nt)];

Phi_rv_tf = [sin(nt)/n, 2/n * (1 - cos(nt)), 0;
          2/n*(cos(nt)-1), 4/n*sin(nt)-3*tf,0;
          0,0,1/n*sin(nt)];

Phi_vr_tf = [3*n*sin(nt),0,0;
          6*n*(cos(nt)-1),0,0;
          0,0,-n*sin(nt)];
Phi_vv_tf = [cos(nt), 2*sin(nt),0
         -2*sin(nt), 4*cos(nt)-3, 0;
          0,0,cos(nt)];

delta_v_0_plus = -inv(Phi_rv_tf) * Phi_rr_tf * delta_r_0;
delta_v_f_minus = Phi_vr_tf * delta_r_0 + Phi_vv_tf * delta_v_0_plus;

% Delta-v at the beginning of the rendezvous maneuver
Delta_v_0 = delta_v_0_plus - delta_v_0_minus;
% Delta-v at the conclusion of the maneuver
Delta_v_f = [0;0;0] - delta_v_f_minus;



% Drawing Reference Chaser Trajector wrt Target
t_vector = linspace(0,tf);
delta_r_t_mat = zeros(length(t_vector),3);
delta_v_t_mat = zeros(length(t_vector),3);

for timestep = 1:length(t_vector)

    t = t_vector(timestep);

    nt = n*t;

    Phi_rr_t = [4 - 3*cos(nt),   0, 0;
              6*(sin(nt) - nt),1,0;
              0,0, cos(nt)];
    
    Phi_rv_t = [sin(nt)/n, 2/n * (1 - cos(nt)), 0;
              2/n*(cos(nt)-1), 4/n*sin(nt)-3*t,0;
              0,0,1/n*sin(nt)];
    
    Phi_vr_t = [3*n*sin(nt),0,0;
              6*n*(cos(nt)-1),0,0;
              0,0,-n*sin(nt)];
    Phi_vv_t = [cos(nt), 2*sin(nt),0
             -2*sin(nt), 4*cos(nt)-3, 0;
              0,0,cos(nt)];

    delta_r_t = Phi_rr_t * delta_r_0 + Phi_rv_t * delta_v_0_plus;
    delta_v_t = Phi_vr_t * delta_r_0 +  Phi_vv_t * delta_v_0_plus;

    delta_r_t_mat(timestep,:) = delta_r_t';
    delta_v_t_mat(timestep,:) = delta_v_t';
end

v_plus_ECI = v_0 + cross(Omega_target, delta_r) + inv(Q_LVLH_ECI)*delta_v_0_plus;

end