function A = project_A(A, params)
%This function projects the input matrix onto the set of interest of
%matrices with 0 on the right bottom K by K square

if isfield(params, 'proj_A')
    switch params.proj_A
        case 'positive'
            A(A<0) = 0;
        case 'negative'
            A(A>0) = 0;
    end
end

A(eye(size(A))>0) = 1;

end