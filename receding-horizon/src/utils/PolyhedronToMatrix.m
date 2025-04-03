function [F, G, b] = PolyhedronToMatrix(X, U)
    % Constraints set X and U in Polyhedron form -> F, G matrix form
    % nc; number of contraint forced by Xc and Uc
    % F; G; % constraints for state and input: Fx+Gu<=1, where 1 is a vector
    Fx = X.A; % Inequalitieis for X only
    Gu = U.A; % Inequalities for U only 
    if numel(Fx)==0
        Fx = zeros(0, X.Dim);
    end
    if numel(Gu)==0
        Gu = zeros(0, U.Dim);
    end
    F = [Fx; zeros(size(Gu, 1), X.Dim)];
    G = [zeros(size(Fx, 1), U.Dim); Gu];
    b = [X.b; U.b]; % Number of Constraints == size(G,1)
end
