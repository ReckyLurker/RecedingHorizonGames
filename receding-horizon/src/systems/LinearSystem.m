classdef LinearSystem < handle
    % Linear Time-Invariant System (LTI) with LQR controller and closed-loop dynamics. 
    % Can be used to express nominal LTI system for robust MPC usecases.
    % utilizes matlab control toolbox 
    
    properties (SetAccess = private)
        % System dynamics: x[k+1] = A*x[k] + B*u[k]
        A;   % State transition matrix
        B;   % Input matrix
        nx;  % Dimension of the state space
        nu;  % Dimension of the input space
        
        % LQR-related properties
        Q;   % State cost matrix for LQR
        R;   % Input cost matrix for LQR
        K;   % LQR feedback matrix: u = K*x
        P;   % Solution to the discrete-time Riccati equation, used in LQR cost function: x'*P*x
        Ak;  % Closed-loop dynamics: A + B*K
    end

    methods (Access = public)

        function obj = LinearSystem(A, B, Q, R)
            % Constructor: Initializes system dynamics and LQR feedback
            % Store system matrices
            obj.A = A;  % State transition matrix
            obj.B = B;  % Input matrix
            obj.Q = Q;  % State cost matrix
            obj.R = R;  % Input cost matrix
            obj.nx = size(A, 1);  % Number of states
            obj.nu = size(B, 2);  % Number of inputs
            % Compute the LQR gain matrix and the cost matrix P
            [K_tmp, obj.P] = dlqr(obj.A, obj.B, obj.Q, obj.R); % from matlab toolbox 
            obj.K = -K_tmp;  % For convinience 
            obj.Ak = obj.A + obj.B * obj.K;  % Closed-loop system matrix
        end

        function x_next = next(obj, x, u)
            % Propagate system dynamics: x[k+1] = A*x[k] + B*u[k]
            x_next = obj.A * x + obj.B * u;  
        end

        function [Ah, Bh] = compute_dynamics(obj, H)
            % Computes the Equality Constraint for the dynamics (X = Ah * x0 + Bh * U)
            Ah = zeros(obj.nx * H, obj.nx);
            Bh = zeros(obj.nx * H, obj.nu);
            for i = 1:H 
                Ah((i-1) * obj.nx + (1:obj.nx), 1:obj.nx) = obj.A ^ i;
                for j = 1:i 
                    Bh((i-1) * obj.nx + (1:obj.nx), (j-1) * obj.nu + (1:obj.nu)) = obj.A ^ (i - j) * obj.B; 
                end
            end
        end

        function Xmpi = compute_MPIset(obj, F, G, b)
            % Computes the Maximal Positively Invariant (MPI) set for the system under LQR control 
            % Function to compute Fpi at time step i
            Fpi = @(i) (F + G * obj.K) * obj.Ak^i;
            % Function to create polyhedron Xpi from Fpi
            Xpi = @(i) Polyhedron(Fpi(i), b);
            % Initialize the MPI set as Xpi(0)
            Xmpi = Xpi(0);
            i = 0;
            % Iteratively compute the MPI set
            while true
                i = i + 1;
                Xmpi_tmp = and(Xmpi, Xpi(i));  % Compute intersection of Xmpi and Xpi(i)
                % Check for convergence (i.e., if the set stops changing)
                if Xmpi_tmp == Xmpi
                    break;
                else
                    Xmpi = Xmpi_tmp;
                end
            end

            if isEmptySet(Xmpi)
                error('MPI Set is Empty');
            else 
                fprintf('MPI Set computed with nu: %d\n', i);
            end
        end
    end
end
