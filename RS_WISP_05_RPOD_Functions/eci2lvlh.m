function delta_r_LVLH = eci2lvlh(delta_r_eci, r_eci, v_eci)
  
    % Reference Frame of the Target
    i_hat = r_eci/norm(r_eci);
    k_hat = cross(r_eci, v_eci) / norm(cross(r_eci, v_eci));
    j_hat = cross(k_hat,i_hat);
    
    % Transformation Matrix from ECI to Space Station Frame
    Q_LVLH_ECI = [i_hat';j_hat';k_hat'];
    delta_r_LVLH = Q_LVLH_ECI * delta_r_eci;
end