function [ us, Us, obj ] = no_need_for_lagrange( Hi, A, k, V, params )

 if strcmp(params.relaxation, 'graph-cut')
     
     [ us, Us, obj ] = lagrange_inner_loop( Hi, A, k, V, params );
     
 else
     
     %Building the quadratic form
     B = [-A, Hi/2 ; Hi'/2, 0];
     D = [eye(V), zeros(V, 1) ; zeros(1, V), 0];
     
     %Linear constraint matrix
     N              = [ones(V, 1), zeros(V, 1) ; 0, 1];
     n_constraints  = size(N, 2);
     [Q, R]         = qr(N);
     
     QBQ      = Q' * B * Q;
     C      = QBQ((n_constraints + 1):end, (n_constraints + 1) :end);
     Gamma  = 2 * (QBQ((n_constraints + 1):end, 1:n_constraints));
     rot    = R(1:n_constraints, 1:n_constraints);
     U1     = (rot \ [2 * k - V , 1]')';
     
     %Need to look where the 0 is after changing the basis (necessarily on one of the two first coordinates)
     Drot   = Q' * D * Q;
     
     %Part corresponding to U1 and U2 (with or without U) respectively
     D1          = Drot(1:n_constraints, 1:n_constraints);
     D2          = Drot((n_constraints + 1):end, (n_constraints + 1):end);
     S           = V - U1 * D1 * U1';
     
     %Linear part of the function
     b           = Gamma * (rot \ [2 * k - V , 1]'); 
     
     
     [U2, ~, obj] = general_spectral_solver(b, -C, S);
     
     U = [U1, U2];
     
     us  = Q*U';
     us  = us(1:end-1)';
     obj = us * Hi - us * A * us';
     Us  = us' * us;
     
 end