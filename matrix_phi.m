% The value G can be written as G = [lambda1 lambda1^2  lambda1^3 ... lambda1^2td-2]  * 
function [d, psi,phi,r] = matrix_phi(t,k)
    d = floor( (2 * k - 1) / t );
    syms delta [1,2 * k - 2]
    syms lambda [1,k]
    syms P
    % generate the matrix AA 
    for i = 1 : k
        % j = 1
        if i ~= 1
            psi(2 * i - 1, 1) = delta(2 * i - 3);
            psi(2 * i, 1) = delta(2 * i - 2);
        end
        % j = 2
        psi(2 * i - 1, 2) = 1 + lambda(i) * P;
        psi(2 * i, 2) = lambda(i) * P;
        % j = 3 to j = 2 * k
        for j = 3 : t * d + 1
            psi(2 * i - 1, j) = lambda(i)^(j - 1);
            psi(2 * i , j) = (j - 1) * lambda(i)^(j - 1);
        end
    end
    psi = psi(1 : t * d + 1 , :);
    % psi is the matrix \Psi = [\Delta , \Lambda_1+P *\Lambda_2, \Lambda_3, ... , \Lambda_{td + 1} \ ]
    % \Psi satisfy that G = -|\Psi|/|\Lambda|, where the matrix \Lambda is the matrix obtained by picking the first td+1 rows of the interpolation matrix.
    
    x = coeffs( det(psi), lambda(1) );
    for i = 1 : 2 * t * d - 2
        for j  = 1 : t * d - 1
            temp = coeffs(x(i), delta(j));
            phi(i,j) = temp(2);
        end
    end
    
    r = rank(phi); % r is the rank of the matrix phi. r should be (t * d - 1) if phi is of full rank. 
end
