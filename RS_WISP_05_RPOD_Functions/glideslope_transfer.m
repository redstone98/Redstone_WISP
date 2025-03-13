function [T, time_vector, Delta_t, Delta_v_m_LVLH_mat, Delta_r_m_LVLH_mat ,delta_r_t_mat, delta_v_t_mat, rho_t_mat, rho_dot_t_mat, delta_v_m_plus_mat] = glideslope_transfer(r_target, v_target, r_chaser, v_chaser, rho_vec_T, rho_dot_0, rho_dot_T, N)

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

% Relative velocity just before launch into the rendezvous trajectory
delta_v_0_minus = Q_LVLH_ECI * delta_v;

% 6. Get rho_0_vec, rho_0, and u_rho

rho_0_vec = delta_r_0 - rho_vec_T;
rho_0 = norm(rho_0_vec);
u_rho = rho_0_vec / rho_0;

% 7. Caculate slope a, total time T, and rho_vec_t, rho_dot_t;

a = (rho_dot_0 - rho_dot_T)/rho_0;
T = 1/a * log(rho_dot_T/rho_dot_0);
Delta_t = T/N;


% 7.1 Generate dataset of rho_t and rho_dot_t
t_vec = linspace(0,T);
rho_t_mat = zeros(length(t_vec),3);
rho_dot_t_mat = zeros(length(t_vec),3);


    for time_index = 1:length(t_vec)
        t = t_vec(time_index);
        rho_t = rho_0 * exp(a*t) + (rho_dot_T/a) * (exp(a*t)-1);
        rho_vec_t = rho_t * u_rho;
        rho_dot_vec_t = (a * rho_t + rho_dot_T) * u_rho;
        rho_t_mat(time_index,:) = rho_vec_t';
        rho_dot_t_mat(time_index,:) = rho_dot_vec_t'; 
    end


% 8.1 Caculate delta_v_m_plus, delta_v_m_minus, Delta_v_m

delta_v_m_minus_mat = zeros(N+1,3);
delta_v_m_minus_mat(1,:) = delta_v_0_minus';
delta_r_m_mat = zeros(N+1,3);
delta_r_m_mat(1,:) = delta_r_0';
delta_v_m_plus_mat = zeros(N,3);
Delta_v_m_LVLH_mat = zeros(N,3);
Delta_r_m_LVLH_mat = zeros(N,3);

delta_r_t_mat = [];
delta_v_t_mat = [];
time_vector = [];

delta_t_vec = linspace(0,Delta_t);
delta_r_t_mat_temp = zeros(length(delta_t_vec),3);
delta_v_t_mat_temp = zeros(length(delta_t_vec),3);


[Phi_rr_Delta_t, Phi_rv_Delta_t, Phi_vr_Delta_t, Phi_vv_Delta_t] = cw_matrix_generator(n,Delta_t);

    for m = 1:N
    
        t_m = m * Delta_t;
        rho_m = rho_0 * exp(a * t_m) + (rho_dot_T / a) * (exp(a*t_m)-1);
        delta_r_m = delta_r_m_mat(m,:)';
        delta_r_m_plus_1 = rho_vec_T + rho_m * u_rho;
        delta_r_m_mat(m+1,:) = delta_r_m_plus_1';
        
        delta_v_m_plus = inv(Phi_rv_Delta_t) * (delta_r_m_plus_1 - Phi_rr_Delta_t * delta_r_m);
        delta_v_m_plus_1_minus = Phi_vr_Delta_t * delta_r_m + Phi_vv_Delta_t * delta_v_m_plus;
        
        delta_v_m_minus_mat(m+1,:) = delta_v_m_plus_1_minus';
        delta_v_m_plus_mat(m,:) = delta_v_m_plus;
        Delta_v_m_LVLH_mat(m,:) =  (delta_v_m_plus - delta_v_m_minus_mat(m,:)')';
        Delta_r_m_LVLH_mat(m,:) = delta_r_m_plus_1';

        for t_index = 1:length(delta_t_vec)
    
            t = delta_t_vec(t_index);
            [Phi_rr, Phi_rv, Phi_vr, Phi_vv] = cw_matrix_generator(n,t);
    
            delta_r_t = Phi_rr * delta_r_m + Phi_rv * delta_v_m_plus;
            delta_v_t = Phi_vr * delta_r_m + Phi_vv * delta_v_m_plus;
    
            delta_r_t_mat_temp(t_index,:) = delta_r_t';
            delta_v_t_mat_temp(t_index,:) = delta_v_t';
        end
    
        delta_r_t_mat = [delta_r_t_mat;delta_r_t_mat_temp];
        delta_v_t_mat = [delta_v_t_mat;delta_v_t_mat_temp];

        time_info = delta_t_vec + (m-1) * Delta_t;

        time_vector = [time_vector;time_info'];
    end

end