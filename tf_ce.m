%% Slightly modified code from Gregory Plett:
% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function [tf_ce, res0, sampling_f, T_len] = tf_ce(SOC_neg, SOC_pos, cse_neg, cse_pos, calculate_ocv_derivative_neg, calculate_ocv_derivative_pos, locs, const)
    M = 10; % Must be at least as big as the statespace.
    sampling_f = 1000;
    T = 1 / sampling_f;
    T_len = 10;
    num_samples = 2 ^ (ceil(log2(sampling_f * T_len)));
    f_vector = 0 : num_samples - 1;
    s = zeros(1, size(f_vector, 2));
    for i = 1 : size(f_vector, 2)
        s(1, i) = 2j * sampling_f * tan(pi * f_vector(i) / num_samples);
    end

    % Copy some variables from cellData structure into more convenient local
    % variable definitions
    L_neg = const.L_neg;
    L_sep = const.L_sep;
    L_pos = const.L_pos;
    L_tot = L_neg + L_sep  + L_pos;
    L_negsep = L_neg + L_sep;
    
    % Check parameters passed to the function
    if(min(locs) < 0 || max(locs) > L_tot+eps)
      error('ERROR (tf_ce): "x" must be in range 0 to %2.2e',L_tot);
    end
    
    % Set up constants and calculations needed by transfer function
    F = const.F;                 % Faraday constant, [Coulomb/mol]
    R = const.R;                 % Gas constant, [J/mol-K]
    T = const.temp;              % Cell temperature, [K]
    t_plus = const.transfer_sep; % Transference number
    zeta = (1 - t_plus) / F;

    Rs_neg = const.radius_neg;       % Particle radius, [m]
    Rs_pos = const.radius_pos;       % Particle radius, [m]
    Ds_neg = const.diffusivity_neg;       % Solid diffusivity [m^2/s]
    Ds_pos = const.diffusivity_pos;       % Solid diffusivity [m^2/s]
    Acell = const.A;   % Current-collector area [m^2]
    as_neg = 3*const.porosity_solid_neg/Rs_neg; % Specific interfacial surf. area
    as_pos = 3*const.porosity_solid_pos/Rs_pos; % Specific interfacial surf. area
    eps1 = const.porosity_electrolyte_neg; % Porosity of negative electrode
    eps2 = const.porosity_electrolyte_sep;
    eps3 = const.porosity_electrolyte_pos; % Porosity of positive electrode
    D1 = const.diffusivity_neg * eps1^const.brug_neg; % Effective ...
    D2 = const.diffusivity_neg * eps2^const.brug_sep; % diffusivities ...
    D3 = const.diffusivity_neg * eps3^const.brug_pos; % of cell regions

    % % Effective conductivities of electrolyte and solid
    kappa_eff_neg = const.kappa_neg;
    kappa_eff_pos = const.kappa_pos;
    sigma_eff_neg = const.sigma_neg;
    sigma_eff_pos = const.kappa_pos;
    
    % % Electrode state of charge for each electrode at this cell SOC               
    theta_neg = SOC_neg;
    theta_pos = SOC_pos;

    % % Copy some variables from cellData structure into more convenient local
    % % variable definitions
    csmax_neg = const.solid_max_c_neg;
    csmax_pos = const.solid_max_c_pos;
    alpha_neg = const.alpha_neg;
    alpha_pos = const.alpha_pos;
    ce0 = const.ce0;
    cs0_neg = csmax_neg * theta_neg; % initial cs in negative
    cs0_pos = csmax_pos * theta_pos; % initial cs in positive
    k_neg = const.kappa_norm_neg/csmax_neg/ce0^(1-alpha_neg); % reaction rate
    k_pos = const.kappa_norm_pos/csmax_pos/ce0^(1-alpha_pos); % ... constants
    % Exchange flux densities at this SOC setpoint
    j0_neg = k_neg * (ce0 * (csmax_neg - cs0_neg))^(1 - alpha_neg) ...
                    * cs0_neg^(alpha_neg);
    j0_pos = k_pos * (ce0 * (csmax_pos - cs0_pos))^(1 - alpha_pos) ...
                    * cs0_pos^(alpha_pos);
    Rct_neg = R*T/(j0_neg*F^2); % Charge-transfer resistance 
    Rct_pos = R*T/(j0_pos*F^2); 
    Rfilm_neg = const.R_film_neg; % Film resistance at particle surface
    Rfilm_pos = const.R_film_pos; 
    Rtot_neg = Rct_neg + Rfilm_neg; % Total resistance across interface
    Rtot_pos = Rct_pos + Rfilm_pos; 

    % % Compute (scaled) transfer function of C_{s,e}(s)/J(s) [Jacobsen-West]
    beta_neg = Rs_neg*sqrt(s./Ds_neg); 
    beta_pos = Rs_pos*sqrt(s./Ds_pos);
    cse_j_neg = Rs_neg/Ds_neg./(1-beta_neg.*coth(beta_neg));
    cse_j_pos = Rs_pos/Ds_pos./(1-beta_pos.*coth(beta_pos));
    % Evaluate derivative of Uocp of both electrodes at given cell SOC
    dudc_neg = calculate_ocv_derivative_neg(cse_neg, const);
    dudc_pos = calculate_ocv_derivative_pos(cse_pos, const);

    % Compute unitless impedance ratio terms for both electrodes
    nu_neg = L_neg*sqrt((as_neg/sigma_eff_neg + as_neg/kappa_eff_neg)...
                     ./(Rtot_neg + dudc_neg.*cse_j_neg/F));
    nu_pos = L_pos*sqrt((as_pos/sigma_eff_pos + as_pos/kappa_eff_pos)...
                     ./(Rtot_pos + dudc_pos.*cse_j_pos/F));

    len = length(locs); locs = locs(:); s = s(:);
    tf_ce = zeros(len,length(s)); % Initialize output to zeros
    % Find eigenvalues that work for separation-of-variables solution
    ce_roots = findLambda(M+1);
    % For each eigenvalue, compute contribution to transfer function
    for n = 1:M,
        lambda = ce_roots(n+1);
        n1 = sqrt(lambda*eps1/D1); Wnn = L_neg*n1; 
        n2 = sqrt(lambda*eps2/D2); Lns = L_neg*n2; Lms = L_negsep*n2;
        n3 = sqrt(lambda*eps3/D3); Wns = L_negsep*n3; Wnt = L_tot*n3;
             Wnp = L_pos*n3;
        
        % Compute scaled coefficients of eigenfunction terms (k1 = 1 for now)     
        k3s = cos(Wnn).*cos(Lns) + D1*n1.*sin(Wnn).*sin(Lns)./(D2*n2);
        k4s = cos(Wnn).*sin(Lns) - D1*n1.*cos(Lns).*sin(Wnn)./(D2*n2);
        k5s = k3s.*(cos(Lms).*cos(Wns)+D2*n2.*sin(Lms).*sin(Wns)./(D3*n3))...
              +k4s.*(sin(Lms).*cos(Wns)-D2*n2.*cos(Lms).*sin(Wns)./(D3*n3));
        k6s = k3s.*(cos(Lms).*sin(Wns)-D2*n2.*sin(Lms).*cos(Wns)./(D3*n3))...
              +k4s.*(sin(Lms).*sin(Wns)+D2*n2.*cos(Lms).*cos(Wns)./(D3*n3));
        % This is used to calculate the value for k1. The 3 terms are the
        % integral of psi^2*epsx over each region. Solved with Mathematica.
        t1 = eps1*(2*Wnn+sin(2*Wnn))/(4*n1);
        t2 = eps2/(4*n2)*(2*(k3s^2+k4s^2)*L_sep*n2+2*k3s*k4s*cos(2*Lns)-...
             2*k3s*k4s*cos(2*Lms)-(k3s-k4s)*(k3s+k4s)*(sin(2*Lns)-sin(2*Lms))); 
        t3 = eps3/(4*n3)*(2*(k5s^2+k6s^2)*L_pos*n3+2*k5s*k6s*cos(2*Wns)-...
             2*k5s*k6s*cos(2*Wnt)-(k5s-k6s)*(k5s+k6s)*(sin(2*Wns)-sin(2*Wnt)));
        k1 = 1/sqrt(t1+t2+t3); % This is proper scaling factor k1
        % Now, properly scale k3, k4, k5, and k6 from k3s (etc).
        k3 = k1*k3s; k4 = k1*k4s; k5 = k1*k5s; k6 = k1*k6s; 
        
        % Now, calculate flux contribution j_n for this eigenvalue
        % Negative electrode
        jn_neg = k1*zeta*nu_neg.*(Wnn*(kappa_eff_neg + ...
             sigma_eff_neg*cosh(nu_neg))*sin(Wnn)+(kappa_eff_neg + ...
             sigma_eff_neg*cos(Wnn)).*sinh(nu_neg).*nu_neg) ./ ...
             (Acell*(kappa_eff_neg+sigma_eff_neg)*(Wnn^2 + ...
             nu_neg.^2).*sinh(nu_neg)); % Value in general
        tf0 = k1*zeta*sin(Wnn)/(Acell*Wnn); % Value when s==0
        jn_neg(s==0) = tf0(:,ones(size(find(s==0))));

        % Positive electrode
        jn_pos = -zeta*nu_pos./(Acell*(kappa_eff_pos+...
             sigma_eff_pos)*(Wnp^2 + nu_pos.^2).*sinh(nu_pos)).*(...
              -k6*Wnp*cos(Wnt)*(sigma_eff_pos+kappa_eff_pos.*cosh(nu_pos)) ...
              +Wnp*(kappa_eff_pos + sigma_eff_pos*cosh(nu_pos)).*...
              (k6*cos(Wns) - k5*sin(Wns)) ...
              +k5*Wnp*(sigma_eff_pos+ ...
              kappa_eff_pos*cosh(nu_pos))*sin(Wnt) ...
              + sinh(nu_pos).*(k5*sigma_eff_pos*cos(Wns)+ ...
              k5*kappa_eff_pos*cos(Wnt) ...
              +k6*sigma_eff_pos*sin(Wns) + ...
              k6*kappa_eff_pos*sin(Wnt)).*nu_pos); % Value in general
        tf0 = -zeta*(k6*(cos(Wns) - cos(Wnt)) + ... % Value when s==0
                           k5*(sin(Wnt)-sin(Wns)))/(Acell*Wnp);
        jn_pos(s==0) = tf0(:,ones(size(find(s==0))));      
        % Compute C_{e,n}(locs,s)/I_app(s) for eigenvalue n
        ce_n = (jn_neg(:) + jn_pos(:))./(s+ce_roots(n+1)); 

        for theLoc = 1:len, % Now, weight by eigenfunction
            x = locs(theLoc); 
            if x < L_neg+eps, % Compute eigenfunction at this location "x"
              Psi = k1*cos(n1*x); % negative electrode region
            elseif x > L_negsep - eps,
              Psi = k5*cos(n3*x)+k6*sin(n3*x); % positive electrode region
            else
              Psi = k3*cos(n2*x)+k4*sin(n2*x); % separator
            end
            % Add the contribution of this eigenvalue to what has already been
            % computed
            tf_ce(theLoc,:) = tf_ce(theLoc,:) + Psi*ce_n.';
            for freq = 1 : size(tf_ce, 2)
                if isnan(tf_ce(1, freq))
                    tf_ce(1, freq) = tf0;
                end
            end
            % Make sure output data structure has "D" term information and the
            % name of the transfer function 
        end
    end
    res0 = zeros(len, 1); % No integrator residue for this transfer function

    

  % Search for "num_roots" roots of Psiprime == 0
  % This is not a very efficient procedure -- could be made a lot faster
  % Requires optimization toolbox for "fzero".  Could replace fzero with
  % a hand-coded bisection algorithm, if this becomes a problem.
  function roots = findLambda(num_roots) 
    roots = 0; w = 2; dL=0.00001; 
    if(num_roots > 1) %if only 1 root then no need to find others!!
      while 1,
        if lambdaFn((w-1)*dL)*lambdaFn(w*dL)<0, % sign change! 
         roots = [roots,fzero(@lambdaFn,[(w-1)*dL w*dL])]; %#ok<AGROW>
         if length(roots)>=num_roots, break; end
        end
        w = w+1; 
      end
    end
  end
   
  % Calculates the derivative of the eigenfunction Psi at the positive-
  % electrode current collector for a proposed eigenvalue lambda. This
  % *should* be zero, and this function is used to search for a zero.
  function Psiprime = lambdaFn(lambda)
    k1 = 1; % True value of k1 is not important when looking for zero
    sle1 = sqrt(lambda*eps1/D1); % shorten names of commonly used variables
    sle2 = sqrt(lambda*eps2/D2); 
    sle3 = sqrt(lambda*eps3/D3);
    k3 = k1*(cos(sle1*L_neg).*cos(sle2*L_neg) + ...
              D1*sle1.*sin(sle1*L_neg).*sin(sle2*L_neg)./(D2*sle2));
    k4 = k1*(cos(sle1*L_neg).*sin(sle2*L_neg) - ...
              D1*sle1.*cos(sle2*L_neg).*sin(sle1*L_neg)./(D2*sle2));

    k5 = k3*(cos(sle2*L_negsep).*cos(sle3*L_negsep) + ...
           D2*sle2.*sin(sle2*L_negsep).*sin(sle3*L_negsep)./(D3*sle3))...
        +k4*(sin(sle2*L_negsep).*cos(sle3*L_negsep) - ...
           D2*sle2.*cos(sle2*L_negsep).*sin(sle3*L_negsep)./(D3*sle3));
    k6 = k3*(cos(sle2*L_negsep).*sin(sle3*L_negsep) - ...
           D2*sle2.*sin(sle2*L_negsep).*cos(sle3*L_negsep)./(D3*sle3))...
        +k4*(sin(sle2*L_negsep).*sin(sle3*L_negsep) + ...
           D2*sle2.*cos(sle2*L_negsep).*cos(sle3*L_negsep)./(D3*sle3));
    Psiprime = -k5.*sle3.*sin(sle3*L_tot) + k6.*sle3.*cos(sle3*L_tot); 
    Psiprime(lambda == 0) = 0;
  end
end
