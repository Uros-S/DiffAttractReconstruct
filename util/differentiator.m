function z = differentiator(signal,t, n_d, n_f, L_0, q, tracking)
% Discrete filtering differentiator with discrete L-adaptation
    
    % Check differentiator parameter validity
    if n_d < 0 || n_f < 0 || n_d + n_f > 12
        disp('Invalid differentiator orders');
        return;
    end

    if tracking && n_d > 4
        disp('Tracking functionality not available for this differentiation order');
        return;
    end

    if ~isnan(q) && tracking
        disp('Discrete L-adaptation not available with tracking');
        return;
    end

    % Utility lambda functions
    sat = @(x) max(-1,min(1,x));
    signmag = @(a,b) abs(a)^b*sign(a);
    
    % Setup differentiator parameters based on n_d and n_f
    switch n_d + n_f
        case 0
            lambda = 1.1;
            k_tau = 2;
        case 1
            lambda = [1.1, 1.5];
            k_tau = 1;
        case 2
            lambda = [1.1, 2.12, 2];
            k_tau = 3;
        case 3
            lambda = [1.1, 3.06, 4.16, 3];
            if n_f == 0
                k_tau = 10;
            else
                k_tau = 6;
            end
        case 4
            lambda = [1.1, 4.57, 9.30, 10.03, 5];
            if n_f == 0
                k_tau = 50;
            else
                k_tau = 10;
            end
        case 5
            lambda = [1.1, 6.75, 20.26, 32.24, 23.72, 7];
            if n_f == 0
                k_tau = 3000;
            else
                k_tau = 400;
            end
        case 6
            lambda = [1.1, 9.91, 43.65, 101.96, 110.08, 47.69, 10];
            if n_f == 0
                k_tau = 1.5e5;
            elseif n_f == 1
                k_tau = 2e4;
            else
                k_tau = 1e4;
            end
        case 7
            lambda = [1.1, 14.13, 88.78, 295.74, 455.40, 281.37, 84.14, 12];
            if n_f == 0
                k_tau = 3e7;
            elseif n_f == 1
                k_tau = 1e7;
            else
                k_tau = 3e6;
            end
        case 8
            lambda = [1.1, 19.66, 171.73, 795.63, 1703.9, 1464.2, 608.99, 120.79, 14];
            if n_f == 0
                k_tau = 1e10;
            elseif n_f == 1
                k_tau = 1.5e10;
            else
                k_tau = 1.5e9;
            end
        case 9
            lambda = [1.1, 26.93, 322.31, 2045.8, 6002.3, 7066.2, 4026.3, 1094.1, 173.72, 17];
            if n_f == 0
                k_tau = 1e12;
            elseif n_f == 1
                k_tau = 7e12;
            elseif n_f == 2
                k_tau = 5e11;
            else
                k_tau = 1.5e11;
            end
        case 10
            lambda = [1.1, 36.34, 586.78, 5025.4, 19895, 31601, 24296, 8908, 1908.5, 251.99, 20];
            if n_f == 0
                k_tau = 1e14;
            elseif n_f == 1
                k_tau = 1e15;
            elseif n_f == 2
                k_tau = 1e13;
            else
                k_tau = 3e12;
            end
        case 11
            lambda = [1.1, 48.86, 1061.1, 12220, 65053, 138954, 143658, 70830, 20406, 3623.1, 386.7, 26];
            if n_f == 0
                k_tau = 8e15;
            elseif n_f == 1
                k_tau = 1e17;
            elseif n_f == 2
                k_tau = 1e15;
            elseif n_f == 3
                k_tau = 8e13;
            else
                k_tau = 2e14;
            end
        case 12
            lambda = [1.1, 65.22, 1890.6, 29064, 206531, 588869, 812652, 534837, 205679, 48747, 6944.8, 623.30, 32];
            if n_f == 0
                k_tau = 3e18;
            elseif n_f == 1
                k_tau = 8e18;
            elseif n_f == 2
                k_tau = 1e16;
            else
                k_tau = 7e15;
            end
    end
    
    % Definition of quasi-continuous rational homogeneous SM controllers
    if tracking && n_d == 0
        alpha = 1.5;
        psi = @(zeta) -sign(zeta(1));
    elseif tracking && n_d == 1
        alpha = 1.5;
        psi = @(zeta) -(signmag(zeta(2),2) + zeta(1)) / ...
                      (zeta(2)^2 + abs(zeta(1)));
    elseif tracking && n_d == 2
        alpha = 4;
        psi = @(zeta) -(zeta(3)^3 + signmag(zeta(2),3/2) + zeta(1)) / ...
                      (abs(zeta(3))^3 + abs(zeta(2))^(3/2) + abs(zeta(1)));
    elseif tracking && n_d == 3
        alpha = 7;
        psi = @(zeta) -(signmag(zeta(4),4) + 2*signmag(zeta(3),2) + 2*signmag(zeta(2),4/3) + zeta(1)) / ...
                      (zeta(4)^4 + 2*zeta(3)^2 + 2*abs(zeta(2))^(4/3) + abs(zeta(1)));
    elseif tracking && n_d == 4
        alpha = 20;
        psi = @(zeta) -(signmag(zeta(5),5) + 6*signmag(zeta(4),5/2) + 5*signmag(zeta(3),5/3) + 3*signmag(zeta(2),5/4) + zeta(1)) / ...
                      (abs(zeta(5))^5 + 6*abs(zeta(4))^(5/2) + 5*abs(zeta(3))^(5/3) + 3*abs(zeta(2))^(5/4) + abs(zeta(1)));
    end
    
    % Initialization result arrays
    z = zeros(length(signal),n_d+1);
    if n_f > 0
        w = zeros(length(signal),n_f);
    end
    if tracking
        zeta = zeros(length(signal),n_d+1);
        zeta(1,:) = 0.0001*ones(1,n_d+1);    % IC cannot be precisely zero since the psi functions would give a NaN
    end
    
    % First iteration setup
    if n_f == 0 && ~tracking
        omega_k = z(1,:)';
    elseif n_f == 0 && tracking
        omega_k = [z(1,:)'; zeta(1,:)'];
    elseif n_f > 0 && ~tracking
        omega_k = [z(1,:)'; w(1,:)'];
    else
        omega_k = [z(1,:)'; w(1,:)'; zeta(1,:)'];
    end    
    
    for k=1:length(signal)-1

        tau = t(k+1)-t(k);
        
        % Initialization iteration variables
        z_next = zeros(n_d+1,1);
        z_k = omega_k(1:n_d+1);
        if n_f > 0
            w_next = zeros(n_f,1);
            if ~tracking
                w_k = omega_k(n_d+1+1:end);
            else
                w_k = omega_k(n_d+1+1:end-(n_d+1));
                zeta_next = zeros(n_d+1,1);
                zeta_k = omega_k(end-n_d:end);
            end
            phi_arg = w_k(1);            % input to D_{n_d,n_f}
        else
            phi_arg = z_k(1)-signal(k); % input to D_{n_d,n_f}
            w_k = z_k(1)-signal(k);     % defined for discrete L-adaptation
        end

        % Discrete L-adaptation
        if isnan(q)
            L = L_0;
            if tracking
                L = L + alpha*L;
            end
        else
            L = L_0*sat(abs(w_k(1))/(L_0*k_tau*tau^(n_d+n_f+1)))^(1+(n_d+n_f+1)*q);
        end

        % D_{n_d,n_f} subsystem
        T_n = zeros(n_d+1,1);

        for i=0:n_d-1
    
            taylor_term = 0;
            
            for j=1:n_d-i
                taylor_term = taylor_term + tau^j/factorial(j)*z_k(j+i+1);
            end

            T_n(i+1) = taylor_term;
            
            if tracking
                zeta_next(i+1) = zeta_k(i+1)+tau*(-lambda(n_d+1-i)*L^((i+1+n_f)/(n_d+n_f+1))*signmag(phi_arg,(n_d-i)/(n_d+n_f+1)) + zeta_k(i+2));
            else
                z_next(i+1) = z_k(i+1)+tau*(-lambda(n_d+1-i)*L^((i+1+n_f)/(n_d+n_f+1))*signmag(phi_arg,(n_d-i)/(n_d+n_f+1))) + taylor_term;
            end
        end
        
        if tracking
            zeta_next(end) = zeta_k(end) + tau*(-lambda(1)*L*sign(phi_arg));
        else
            z_next(end) = z_k(end) + tau*(-lambda(1)*L*sign(phi_arg));
        end
    
        % \Omega_{n_d,n_f} subsystem
        if n_f > 0
            for i=1:n_f-1
                w_next(i) = w_k(i) + tau*(-lambda(end+1-i)*L^(i/(n_d+n_f+1))*signmag(w_k(1),(n_d+n_f+1-i)/(n_d+n_f+1)) + w_k(i+1));
            end
            
            if tracking
                w_next(end) = w_k(end) + tau*(-lambda(n_d+2)*L^(n_f/(n_d+n_f+1))*signmag(w_k(1),(n_d+1)/(n_d+n_f+1)) + zeta_k(1)-z_k(1)+signal(k));
            else
                w_next(end) = w_k(end) + tau*(-lambda(n_d+2)*L^(n_f/(n_d+n_f+1))*signmag(w_k(1),(n_d+1)/(n_d+n_f+1)) + z_k(1)-signal(k));
            end
        
        end

        % Tracking subsystem
        if tracking
            for i=0:n_d-1
                z_next(i+1) = z_k(i+1) + T_n(i+1);
            end
            z_next(end) = z_k(end) + tau*alpha*L_0*psi(1/L_0*zeta_k);
        end
        
        % Save iteration results
        z(k+1,:) = z_next';
        if n_f > 0
            w(k+1,:) = w_next';
        end
        if tracking
            zeta(k+1,:) = zeta_next';
        end
        
        % Setup variables for next iteration
        if k < (length(signal)-1)
            if n_f == 0 && ~tracking
                omega_k = z_next;
            elseif n_f > 0 && ~tracking
                omega_k = [z_next; w_next];
            else
                omega_k = [z_next; w_next; zeta_next];
            end
        end
    end

end