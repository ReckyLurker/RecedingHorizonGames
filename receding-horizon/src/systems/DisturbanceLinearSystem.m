classdef DisturbanceLinearSystem < LinearSystem
    % Linear System with Additive Uncertainty: x_{k+1} = Ax_k + Bu_k + w_k
    % Models a linear system with "only" additive disturbances
    % Utilizes MPT3 toolbox.
    
    properties (SetAccess = private)
        W % Convex set representing the disturbance
        Z % Minimal Robust Positively Invariant (mRPI) set\
        pd; % Distribution of Sampling (Assuming normal distribution N(0.5, 0.5) to sample from 0 to 1)
    end

    methods (Access = public)

        function obj = DisturbanceLinearSystem(A, B, Q, R, D, W)
            % Constructor for DisturbanceLinearSystem
            % A, B, Q, R: System matrices (inherited from LinearSystem)
            % D: Transformation matrix for disturbance
            % W: Convex set of disturbances

            obj = obj@LinearSystem(A, B, Q, R);  % Call superclass constructor
            
            if(size(D, 1) == 0) 
                obj.W = W;
            else 
                obj.W = D * W;
            end
            
            obj.Z = obj.compute_mrpi_set(1e-4);  % Compute the mRPI set upto a given tolerance

            mu = 0.5;   % Mean
            sigma = 1; % Standard deviation
            obj.pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
            obj.pd = truncate(obj.pd, 0, 1);
        end

        function [x_next,disturbance] = next(obj, x, u)
            % Propagate the system dynamics with a random disturbance
            % x: Current state, u: Control input
            % Returns the next state x_new
            disturbance = obj.pick_random_disturbance();  % Pick a random disturbance
            x_next = next@LinearSystem(obj, x, u) + disturbance;  % Add disturbance to next state
        end

    end

    methods (Access = public)

        function disturbance = pick_random_disturbance(obj)
            % Pick a random disturbance from a uniform distribution within the set W
            vertices = obj.W.V;  % Get vertices of the convex set W
            max_disturbance = max(vertices)';  % Maximum bounds of the disturbance
            min_disturbance = min(vertices)';  % Minimum bounds of the disturbance

            % Generate a random disturbance inside the convex set W
            while true
                disturbance = random(obj.pd, obj.nx, 1) .* (max_disturbance - min_disturbance) + min_disturbance;  % Generate random disturbance
                if obj.W.contains(disturbance)  % Check if the disturbance is inside W
                    break
                end
            end
        end

        function mRPI_set = compute_mrpi_set(obj, epsilon)
            % Computes an approximation of the Minimal Robust Positively Invariant (mRPI) set
            % x^{+} = Ax + w with w ∈ W
            % Algorithm-I from Rakovic et al., with tolerance 'epsilon'
            
            [nx, ~] = size(obj.Ak);  % Get the state dimension from matrix A
            iteration_count = 0;  % Initialize iteration counter
            alpha = 1000;  % Initial alpha value for the algorithm
            max_support_sum = 1000;  % Initialize upper bound for the support function sum

            while(alpha > epsilon/(epsilon + max_support_sum))
                iteration_count = iteration_count + 1;
                alpha = max(obj.W.support(obj.Ak^iteration_count * obj.W.A') ./ obj.W.b);  % Compute support function ratio
                support_sums = zeros(2 * nx, 1);  % Initialize support sum vector

                % Sum the support functions over the iterations
                for i = 1:iteration_count
                    support_sums = support_sums + obj.W.support([obj.Ak^i, -obj.Ak^i]);
                end

                max_support_sum = max(support_sums);  % Update the max support sum
            end
            fprintf('The values of r: %d, rho: %.4f\n', iteration_count, alpha);
            % Construct the mRPI set by summing the transformed sets
            mRPI_set = obj.W;
            for i = 1:iteration_count - 1
                mRPI_set = mRPI_set + obj.Ak^i * obj.W;  % Sum the set transformations
            end
            mRPI_set = (1 / (1 - alpha)) * mRPI_set;  % Apply scaling factor based on alpha
        end
    end
end
