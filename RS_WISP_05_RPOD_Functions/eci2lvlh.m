function delta_r_LVLH = eci2lvlh(delta_r_eci, r_eci, v_eci)
  
    % Reference Frame of the Target
    i_hat = r_eci/norm(r_eci);
    j_hat = v_eci/norm(v_eci);
    k_hat = cross(i_hat, j_hat);
    
    % Transformation Matrix from ECI to Space Station Frame
    Q_LVLH_ECI = [i_hat';j_hat';k_hat'];
    delta_r_LVLH = Q_LVLH_ECI * delta_r_eci;
end