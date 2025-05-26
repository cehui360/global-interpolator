function [Dx, Dy, Dxx, Dyy, Dxy] = Kronecker5(N, M, t)
    % Checking input parameters
    if N < 2
        error('N must be greater than or equal to 2');
    end
    if M < 2
        error('M must be greater than or equal to 2');
    end
    if t <= 0
        error('The grid spacing t must be positive');
    end

    % Identity matrix
    I_M = speye(M);
    I_N = speye(N);

    % Construct the first and second derivative matrices
    A1_M = diff_matrix(M, 1); 
    A1_N = diff_matrix(N, 1); 
    A2_M = diff_matrix(M, 2); 
    A2_N = diff_matrix(N, 2); 

    % The Kronecker product generates the derivative matrix
    Dx = (1 / t) * kron(I_M, A1_N);
    Dy = (1 / t) * kron(A1_M, I_N);
    Dxx = (1 / t^2) * kron(I_M, A2_N);
    Dyy = (1 / t^2) * kron(A2_M, I_N);
    Dxy = Dy * Dx;

    function A = diff_matrix(K, order)
        % Building the difference matrix
        if K < 2
            error('The matrix dimension K must be greater than or equal to 2');
        end
        switch order
            case 1 % First derivative
                A = spdiags([-ones(K, 1), ones(K, 1)], [-1, 0], K, K);
            case 2 % Second derivative
                A = spdiags([ones(K, 1), -2 * ones(K, 1), ones(K, 1)], [-1, 0, 1], K, K);
            otherwise
                error('Only first or second derivatives are supported');
        end
        if order == 2 && K > 2
            A(1, 1) = 0; A(1, 2) = 0;
            A(K, K) = 0; A(K, K-1) = 0;
        end
    end
end
