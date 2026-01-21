function re = two_kinetic_model(Csolid, alpha, rho, F, parameters)
    % Langmuir model with parameters depending on control variables (T, P, F)
    % No derived properties (ρ, RE) - directly models in TPF space
    %
    % Model: Linear dependencies of K_m, k_max, and k_diff on (T, P, F)
    %   K_m(T, P) - Half-saturation depends on temperature and pressure
    %   k_max(T, P, F) - Maximum rate depends on all three conditions
    %   k_diff(T, P) - Diffusion rate depends on temperature and pressure
    %
    % This breaks the k_max/k_diff correlation by giving them different
    % functional forms (k_max depends on F, k_diff does not)

    import casadi.*

    %% Extract parameters (indices 44-54, 11 parameters)
    k_w_0      = parameters{44};  % Base half-saturation [kg/m3]
    a_w        = parameters{45};  % K_m temperature coefficient [kg/m3/K]
    b_w        = parameters{46};  % K_m pressure coefficient [kg/m3/bar]
    k_d_0      = parameters{47};  % Base maximum rate [kg/m3/s]
    a_d        = parameters{48};  % k_max temperature coefficient [kg/m3/s/K]
    %k_max_P   = parameters{49};  % k_max pressure coefficient [kg/m3/s/bar]
    %k_max_F   = parameters{50};  % k_max flow coefficient [kg/m3/s/(kg/s)]
    %k_diff_0  = parameters{51};  % Base diffusion rate [1/s]
    %k_diff_T  = parameters{52};  % k_diff temperature coefficient [1/s/K]
    %k_diff_P  = parameters{53};  % k_diff pressure coefficient [1/s/bar]
    %n_order   = parameters{54};  % Extraction order [-]

    
    %%
    k_w = k_w_0 * (rho / 800)^a_w * (F / 5)^b_w * 1e-13;
    k_d = k_d_0 * (rho / 800)^a_d * 1e-13; 

    %%
    re = (k_w * alpha + k_d * (1-alpha) ) * Csolid;

end
