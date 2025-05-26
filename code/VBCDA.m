function [s, z, u] = VBCDA(s0, z0, u0,para,m,testdata,minx,miny,r,c) 
    % Initialize variables
    [N,M] = size(u0);
    g = reshape(u0,N*M,1);
    s = reshape(s0,N*M,1);
    z = reshape(z0,N*M,1);
    u = reshape(u0,N*M,1);
    e = ones(N,M); 
    e = reshape(e,N*M,1);
    k = 0;
    % Parameter setting
    gamma_s = 1;
    gamma_z = 1;
    gamma_u = 1.5;
    xi = 0.001;
    max_iter = para.max_iter;
    alpha = para.alpha; 
    beta = para.beta; 
    delta=para.delta;
    epsilon = para.epsilon;
    tolerance = 10^(-3);
    I = speye(N*M,N*M);
    while k < max_iter
        % Finite difference matrix
        [Dxv,Dyv,Dxxv,Dyyv,Dxyv] = second_order_derivatives(N,M,m,u0);
        
        % Differential operator
        Delta1 = Dxxv.^2 + Dyyv.^2 + 2.*Dxyv.^2;
        Delta1 = reshape(Delta1,N*M,1);
        R_Delta1 = spdiags(Delta1,0,N*M,N*M);
        Delta2 = Dxv.^2 + Dyv.^2;
        Delta2 = reshape(Delta2,N*M,1);
        R_Delta2 = spdiags(Delta2,0,N*M,N*M);

        % Compute Dx Dy Dxx Dyy Dxy 
        [Dx,Dy,Dxx,Dyy,Dxy] = Kronecker5(r,c,m);
        
        % Compute A_s\b_s,A_z\b_z
        A_s = 2*xi*R_Delta2 + 2*epsilon*(alpha - beta)*(Dx'*Dx + Dy'*Dy) + (alpha - beta)/(2*epsilon)*I; % matrix for s
        A_z = 2*delta*R_Delta1 + 2*epsilon*beta*(Dx'*Dx + Dy'*Dy) + beta/(2*epsilon)*I; % matrix for z
        b_s = (alpha - beta)/(2*epsilon)*e; % vector for s
        b_z = beta/(2*epsilon)*e; % vector for z

        % Compute the search directions d_s^k and d_z^k
        d_s_k = compute_search_direction_s(A_s, b_s, s);
        d_z_k = compute_search_direction_z(A_z, b_z, z);

        % Compute alpha_s^k and alpha_z^k
        alpha_s_k = gamma_s * (-A_s * s - b_s)' * d_s_k / (d_s_k' * A_s * d_s_k);
        alpha_z_k = gamma_z * (-A_z * z - b_z)' * d_z_k / (d_z_k' * A_z * d_z_k);
        
        % Update s and z
        s_new = s + alpha_s_k * d_s_k;
        z_new = z + alpha_z_k * d_z_k;

        % Compute Rz Rs
        Rs = spdiags(s_new.^2,0,N*M,N*M);
        Rz = spdiags(z_new.^2,0,N*M,N*M);
        
        % Compute A_s\b_s
        A_u = 2*delta*(Dxx'*Rz*Dxx + Dyy'*Rz*Dyy + 2.*Dxy'*Rz*Dxy) +2*xi*(Dx'*Rs*Dx + Dy'*Rs*Dy) + 2*I;
        b_u = 2*g;

        % Compute the search direction d_u^k_b
        d_u_k = compute_search_direction_u(A_u, b_u, u);

        % Compute alpha_u^k
        alpha_u_k = gamma_u * (-A_u * u - b_u)' * d_u_k / (d_u_k' * A_u * d_u_k);

        % Update u_b
        u_new = u + alpha_u_k * d_u_k;

        k = k + 1;
        % Checking convergence
        disp(['iterations: ' num2str(k)])
        if norm(u_new - u) < tolerance
            break;
        end
        if k > 100
            break
        end
        % Updating variables
        s = s_new;
        z = z_new;
        u = u_new;
        zs = reshape(u,N,M);
    end
end

function d_s_k = compute_search_direction_s(A_s, b_s, s)
    % Implement the computation of the search direction for s
    % Placeholder implementation
    d_s_k = -(A_s*s - b_s); % Example placeholder
end

function d_z_k = compute_search_direction_z(A_z, b_z, z)
    % Implement the computation of the search direction for z
    % Placeholder implementation
    d_z_k = -(A_z*z - b_z); % Example placeholder
end

function d_u_k = compute_search_direction_u(A_u, b_u, u)
    % Implement the computation of the search direction for u
    % Placeholder implementation
    d_u_k = -(A_u*u - b_u); % Example placeholder
end
