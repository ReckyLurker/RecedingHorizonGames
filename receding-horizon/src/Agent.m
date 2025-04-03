classdef Agent < DisturbanceLinearSystem 
    % Defines Agent class for all active agents of the battery charging game.

    properties(SetAccess = private)
        id; % Player id 
        H; % Finite horizon length 

        % Player Polyhedrons [Local]
        X_nominal; 
        U_nominal; 
        Xmpi_nominal; % Terminal Constraint Set 
        Z_shift; 

        % Initial State 
        x_init; 

        % Nominal System Dynamics [over the horizon, precomputed] 
        Ah; 
        Bh; 

        % Cost Function Parameters 
        f;
        mu; 

        % Cost Function [Local Dynamics, Cost Functions, etc., precomputed]
        Qh; % J = X'QhX + Rh'URh + Rh'(U-)Rh + fh'*U 
        Rh; 
        fh; 
        Ph; % Lyapunov Suboptimality 
        suboptimal_lyapunov_cost;

        % Shared Constraints 
        alpha; % Fraction of available global/shared resources   
        F_global;
        G_global;
        b_global;
        b_reduce; 

        % Steady State 
        SS; % Projection Polyhedron
        ss_state; 
        ss_lambda; 
        ss_z; 
        ss_tstate;
        ss_tlambda;
        ss_tz;
        ss_tau;
        ss_nu;
        ss_sigma;
        ss_Ga;
        ss_Gb;


        % Global Inequalities Ax >= b form (because of the operator spilitting algorithm requirements) 
        % Agent knows only it's local block of the global constraint
        Ga;
        Gb; 

        % vGNE Optimizer steps [updated automatically]
        TAU; 
        NU;
        SIGMA; 

        state; % Player State
        lambda; % Lagrange Multiplier 
        z; % Player Consensus Variable

        % Temporary variables 
        tstate; % To store next player state
        tlambda; % To store next lagrange multipliers
        tz; % Player Consensus Variable 

        U_player; % Player State Polyhedron [Used for projection]
    end

    methods(Access = public)
        function obj = Agent(id, H, A, B, Q, R, f, mu, X, U, W, x_init, F_global, G_global, b_global, alpha)
            tic;  
            obj = obj@DisturbanceLinearSystem(A, B, Q, R, [], W);
            
            obj.id = id; 
            obj.H = H; 
            obj.f = f; 
            obj.mu = mu; 
            
            % Shared constraints 
            obj.alpha = alpha; 
            obj.F_global = F_global;
            obj.G_global = G_global;
            obj.b_global = b_global; 

            % Initial State 
            obj.x_init = x_init; 

            % Tighten Local Constraints 
            obj.tighten_local_constraints(X, U);

            % Compute Shared Constraint Reduction 
            obj.compute_shared_constraint_tightening();


            % Formulate Nominal System Dynamics 
            [obj.Ah, obj.Bh] = obj.compute_dynamics(H);

            % Formulate Cost Function 
            [obj.Qh, obj.Rh, obj.fh] = obj.compute_cost_function(f);

            % SDP Lyapunov Suboptimality Constraint 
            obj.construct_SDP_constraint();
            
            time_taken = toc; 
            fprintf('Agent constructed in %.6f sec\n', time_taken);
        end

        function tighten_local_constraints(obj, X, U)
            obj.X_nominal = X - obj.Z; 
            obj.U_nominal = U - obj.K*obj.Z; 
        end

        function compute_shared_constraint_tightening(obj)
            % Tighten b_global 
            e = sdpvar(obj.nx, 1);
            C = [obj.Z.A*e <= obj.Z.b];
            Fk = obj.F_global + obj.G_global * obj.K; 
            options = sdpsettings('solver', 'sedumi', 'verbose', 0); 
            obj.b_reduce = zeros(size(obj.F_global, 1), 1);
            for j = 1:size(obj.F_global, 1)
                J = Fk(j, :) * e;
                optimize(C, -J, options); 
                obj.b_reduce(j) = value(J);
            end
        end

        function delta = communicate_shared_constraint_tightening(obj)
            delta = obj.b_reduce; 
        end

        function update_shared_constraint(obj, delta)
            obj.b_global = obj.b_global - delta; 
        end

        function compute_steady_state_parameters(obj)
            obj.SS = obj.X_nominal * obj.U_nominal; 
            ss_A = [obj.SS.A; (eye(obj.nx)-obj.A) -obj.B; (obj.A-eye(obj.nx)) obj.B]; % Along with inequality constraint 
            ss_b = [obj.SS.b; zeros(obj.nx,1); zeros(obj.nx, 1)];

            obj.SS = Polyhedron(ss_A, ss_b);

            obj.ss_state = zeros(obj.nx + obj.nu, 1);
            obj.ss_lambda = rand(size(obj.F_global, 1), 1); % lagrange multiplier 
            obj.ss_z = zeros(size(obj.b_global, 1), 1); % Consensus Variable 

            obj.ss_Ga = [obj.F_global obj.G_global];
            obj.ss_Gb = [obj.b_global];

            delta = 150; % Hardcoded 
            A1 = 0.0;
            A2 = 0.0;
            for i = 1:size(obj.ss_Ga, 2) 
                A1 = max(A1, sum(abs(obj.ss_Ga(:, i))));
            end
            for i = 1:size(obj.ss_Ga, 1)
                A2 = max(A2, sum(abs(obj.ss_Ga(i, :))));
            end
            obj.ss_tau = 0.75 / (A1 + delta);
            obj.ss_sigma = 0.75 / (A2 + 2 + delta); % Assuming weights are 1/N for each of the players. 
            obj.ss_nu = 0.75 / (2 + delta); % Assuming weights are 1/N for each of the players. 
            fprintf('Choosing Steady State Step Sizes, TAU: %.6e , NU: %.6e, SIGMA: %.6e for player %d\n', obj.ss_tau, obj.ss_nu, obj.ss_sigma, obj.id);
        end

        function update_ss_state(obj, N, players)
            weight = 1/N; % In the graph 
            ss_s = obj.ss_state(1:obj.nx);
            ss_u = obj.ss_state(obj.nx+1:obj.nx+obj.nu);

            pseudo_grad_s = 2 * ss_s' * obj.Q; 
            pseudo_grad_u_local = obj.f' + 2 * ss_u' * obj.R;
            pseudo_grad_u_coupled = zeros(1,obj.nu);
            
            for i = 1:N 
                if i == obj.id 
                    continue; 
                end 
                ss_u_i = players{i}.ss_state(players{i}.nx+1:players{i}.nx+players{i}.nu);
                pseudo_grad_u_coupled = pseudo_grad_u_coupled + 2 * obj.mu * ss_u_i' * obj.R; 
            end

            pseudo_grad = [pseudo_grad_s' ; pseudo_grad_u_local' + pseudo_grad_u_coupled']; 

            obj.ss_tstate = obj.ss_state - obj.ss_tau*(pseudo_grad + obj.ss_Ga' * obj.ss_lambda);

            temp = obj.SS.project(obj.ss_tstate);
            obj.ss_tstate = temp.x; 
            obj.ss_tz = zeros(size(obj.ss_Ga, 1), 1);
            for i = 1:N 
                if i == obj.id 
                    continue; % 0 there 
                end
                [~, player_lambda, ~] = players{i}.communicate_ss_state();
                obj.ss_tz = obj.ss_tz + obj.ss_nu*weight*(obj.ss_lambda - player_lambda);
            end         
        end

        function update_ss_lambda(obj, N, players)
            weight = 1/N; 
            obj.ss_tlambda = obj.ss_lambda + obj.ss_sigma*(obj.ss_Ga * (2*obj.ss_tstate - obj.ss_state) - obj.ss_Gb);
            for i = 1:N 
                if i == obj.id 
                    continue; 
                end
                [~, player_lambda, player_z] = players{i}.communicate_ss_state();
                [~, ~, player_z_update] = players{i}.communicate_ss_update();
                obj.ss_tlambda = obj.ss_tlambda + obj.ss_sigma*weight*(2*(obj.ss_tz - player_z_update) - (obj.ss_z - player_z));
                obj.ss_tlambda = obj.ss_tlambda + obj.ss_sigma*weight*(obj.ss_lambda - player_lambda);
            end
            obj.ss_tlambda = max(obj.ss_tlambda, 0); % Project to R+
        end

        function [ss_state, ss_lambda, ss_z] = communicate_ss_state(obj);
            ss_state = obj.ss_state;
            ss_lambda = obj.ss_lambda;
            ss_z = obj.ss_z; 
        end

        function [ss_tstate, ss_tlambda, ss_tz] = communicate_ss_update(obj)
            ss_tstate = obj.ss_tstate;
            ss_tlambda = obj.ss_tlambda;
            ss_tz = obj.ss_tz; 
        end

        function finish_ss_update(obj)
            obj.ss_state = obj.ss_tstate;
            obj.ss_lambda = obj.ss_tlambda;
            obj.ss_z = obj.ss_tz;
        end

        function shift_system(obj)
            obj.X_nominal = -obj.ss_state(1:obj.nx) + obj.X_nominal; 
            obj.U_nominal = -obj.ss_state(obj.nx+1:obj.nu+obj.nx) + obj.U_nominal; 
            obj.Z_shift = obj.ss_state(1:obj.nx) + obj.Z; 
        end

        function construct_terminal_set(obj)
            % Compute Terminal Constraint Set (A Positively Invariant Set)
            [F, G, b] = PolyhedronToMatrix(obj.X_nominal, obj.U_nominal);
            F = [F; obj.F_global];
            G = [G; obj.G_global];
            b = [b; obj.alpha * obj.b_global];
            obj.Xmpi_nominal = obj.compute_MPIset(F, G, b);  
        end

        function construct_SDP_constraint(obj)
            obj.Ph = [obj.Ah' * obj.Qh * obj.Ah obj.Ah' * obj.Qh * obj.Bh; 
                            obj.Bh' * obj.Qh * obj.Ah obj.Rh + obj.Bh'*obj.Qh*obj.Bh];
        end

        function update_SDP_constraint(obj)
            obj.suboptimal_lyapunov_cost = obj.compute_suboptimal_cost();
        end

        function cost = compute_suboptimal_cost(obj)
            s_H = obj.state(1:obj.nx);
            obj.state(1:obj.nx) = obj.A*obj.state(1:obj.nx) + obj.B*obj.state(obj.nx+1:obj.nx+obj.nu); % s^v_{1|k}
            for i = 0:obj.H-2
                s_H = obj.A*s_H + obj.B*obj.state(obj.nx + 1 + i*obj.nu:obj.nx + (i+1)*obj.nu);
                obj.state(obj.nx + 1 + i*obj.nu:obj.nx + (i+1)*obj.nu) = obj.state(obj.nx + 1 + (i+1)*obj.nu:obj.nx + (i+2)*obj.nu); % u^v_{i|k+1} ← u^v_{i+1|k}
            end
            obj.state(obj.nx + 1 + (obj.H-1)*obj.nu:obj.nx + (obj.H)*obj.nu) = obj.K * s_H; 

            cost = obj.state' * obj.Ph * obj.state; 
        end

        function [x_next, u, disturbance] = next(obj, x)
            % Compute the control input 
            s0 = obj.state(1:obj.nx) + obj.ss_state(1:obj.nx); 
            u_nominal = obj.ss_state(obj.nx+1:obj.nx+obj.nu) + obj.state(obj.nx+1:obj.nx+obj.nu); % state = [s0, u_n0, u_n1, ..., u_n(H-1)]
            u = u_nominal + obj.K * (x - s0);
            [x_next,disturbance] = next@DisturbanceLinearSystem(obj, x, u);
            ss_diff = s0 - obj.ss_state(1:obj.nx);  
            fprintf('SS Norm Player %d: %.6f\n', obj.id, norm(ss_diff) / (norm(obj.ss_state(1:obj.nx) + 1e-6)));
        end

        function update_polyhedron(obj, x_init)
            obj.x_init = x_init; 
            obj.construct_polyhedron();
        end

        function prepare_sim(obj)
            obj.construct_polyhedron();
            obj.construct_global_ineq();
            obj.compute_step_sizes();
            obj.state = zeros(obj.nx + obj.H*obj.nu, 1);
            obj.lambda = rand(size(obj.Ga, 1), 1); % lagrange multiplier 
            obj.z = zeros(size(obj.Ga, 1), 1); % Consensus Variable 
        end

        % vGNE Solution Algorithms [noncooperative in cost reduction, coopertative in constraint satisfaction]
        % Algorithm adapted from Algorithm-I of https://doi.org/10.1016/j.automatica.2019.01.008
        function update_state(obj, N, players, t)
            weight = 1/N; 
            pseudo_grad = obj.oracle(N, players);
            obj.tstate =  obj.state - obj.TAU*(pseudo_grad + obj.Ga'*obj.lambda);
            if t == 1 
                temp = obj.U_player.project(obj.tstate); % Without Dissipativity 
                obj.tstate = temp.x; 
            else
                obj.tstate = obj.project_point(obj.tstate); % With Dissipativity
            end
            obj.tz = zeros(size(obj.Ga,1), 1);
            for i = 1:N 
                if i == obj.id 
                    continue; % 0 there 
                end
                [~, player_lambda, ~] = players{i}.communicate_state();
                obj.tz = obj.tz + obj.NU*weight*(obj.lambda - player_lambda);
            end
        end

        function x_proj = project_point(obj, x0)
            % tic; 
            x = sdpvar(obj.nx+obj.H*obj.nu,1);
            objective = (x-x0)' * (x-x0);
            constraints = [x' * obj.Ph * x <= obj.suboptimal_lyapunov_cost];
            constraints = [constraints; obj.U_player.A * x <= obj.U_player.b];
            options = sdpsettings('verbose', 0, 'solver', 'mosek');
            optimize(constraints, objective, options);
            x_proj = value(x);
            % fprintf('Time Taken for Optimization: %.6f\n', toc);
            % if ~any(isnan(x_proj))
            %     disp('Projection successful');
            % else
            %     warning('Projection failed - no feasible solution found');
            %     x_proj = NaN(n,1);
            % end
        end

        function update_lambda(obj, N, players)
            weight = 1/N; 
            obj.tlambda = obj.lambda + obj.SIGMA*(obj.Ga * (2*obj.tstate - obj.state) - obj.Gb);
            for i = 1:N 
                if i == obj.id 
                    continue; 
                end
                [~, player_lambda, player_z] = players{i}.communicate_state();
                [~, ~, player_z_update] = players{i}.communicate_update();
                obj.tlambda = obj.tlambda + obj.SIGMA*weight*(2*(obj.tz - player_z_update) - (obj.z - player_z));
                obj.tlambda = obj.tlambda + obj.SIGMA*weight*(obj.lambda - player_lambda);
            end
            obj.tlambda = max(obj.tlambda, 0); % Project to R+
        end

        function finish_update(obj)
            obj.state = obj.tstate;
            obj.lambda = obj.tlambda;
            obj.z = obj.tz;
        end
        
        % To communicate State 
        function [state, lambda, z] = communicate_state(obj)
            state = obj.state; 
            lambda = obj.lambda; 
            z = obj.z; 
        end

        function [tstate, tlambda, tz] = communicate_update(obj)
            tstate = obj.tstate;
            tlambda = obj.tlambda;
            tz = obj.tz;
        end

        function [J,J_stage] = predict_cost(obj, N, players)
            % Nominal System Cost (which was used for designing LQR Controller)
            s0 = obj.state(1:obj.nx);
            u0 = obj.state(obj.nx+1:obj.nx+obj.nu);
            u_nominal_seq = obj.state(obj.nx+1:obj.nx+obj.H*obj.nu);

            J_nominal = s0' * (obj.Q + obj.Ah' * obj.Qh * obj.Ah) * s0 + 2 * s0' * obj.Ah' * obj.Qh * obj.Bh * u_nominal_seq + ...
                        + u_nominal_seq' * (obj.Rh + obj.Bh' * obj.Qh * obj.Bh) * u_nominal_seq; 
            J_linear = obj.fh' * u_nominal_seq;
            
            J_coupled = 0.0;

            J_nominal_stage = s0' * obj.Q * s0 + u0' * obj.R * u0;
            J_linear_stage = obj.f' * u0; 
            J_coupled_stage = 0; 
            for i = 1:N 
                if i == obj.id 
                    continue; % Already Considered in Nominal Cost 
                end
                [player_state, ~, ~] = players{i}.communicate_state();
                u_i_nominal_seq = player_state(obj.nx + 1:obj.nx + obj.H*obj.nu);
                J_coupled = J_coupled + 2*(u_i_nominal_seq)' * obj.Rh * (u_nominal_seq);
                J_coupled_stage = J_coupled_stage + obj.mu * u0' * obj.R * u_i_nominal_seq(1:obj.nu);
            end
            J = J_nominal + J_linear + J_coupled;
            J_stage = J_nominal_stage + J_linear_stage + J_coupled_stage; 
        end
    end

    methods(Access = private)
         function [Qh, Rh, fh] = compute_cost_function(obj, f)
            Qh = [];
            Rh = obj.R;
            fh = f;
            for i = 1:obj.H-1
                Qh = blkdiag(Qh, obj.Q);
                Rh = blkdiag(Rh, obj.R);
                fh = [fh; f];
            end
            Qh = blkdiag(Qh, obj.P);
        end

        function construct_polyhedron(obj)
            % Computes the Polyhedron form for the nominal system with tightened constraints (Rigid Tube Robust MPC)
            % Equality Constraints (Dynamics) are substituted to retain only the inequality constraints
            % Xinit = x_init + obj.Z; 

            % if isEmptySet(Xinit)
            %     error('Initial Condition is outside nominal set');
            % end

            % x in X_nominal and u in U_nominal 
            [F, G, b] = PolyhedronToMatrix(obj.X_nominal, obj.U_nominal);
            nc = size(F, 1);

            % Initial Constraint (x in Xinit)
            Fh0 = [-obj.Z_shift.A; F]; 
            bh = [obj.Z_shift.b - obj.Z_shift.A*obj.x_init; b];
            nc0 = size(obj.Z_shift.A, 1);

            % Over the Horizon
            Fh = [];
            Gh = G; 
            for i = 1:obj.H-1 
                Fh = blkdiag(Fh, F);
                Gh = blkdiag(Gh, G);
                bh = [bh; b];
            end

            % Terminal Constraint (x_h in Xmpi_nominal)
            nch = size(obj.Xmpi_nominal.A, 1); 
            Fh = blkdiag(Fh, obj.Xmpi_nominal.A);
            bh = [bh; obj.Xmpi_nominal.b];

            % Final Assembly
            Fh0 = [Fh0; zeros(obj.H * nc - nc + nch, obj.nx)];
            Fh = [zeros(nc0 + nc, obj.H*obj.nx); Fh];
            Gh = [zeros(nc0, obj.H * obj.nu); Gh; zeros(nch, obj.H*obj.nu)];

            % Construction of Polyhedron
            U_A_x0 = Fh0 + Fh * obj.Ah; 
            U_A_u = Fh * obj.Bh + Gh; 
            U_A = [U_A_x0 U_A_u];
            U_b = bh; 
            obj.U_player = Polyhedron(U_A, U_b);

            if isEmptySet(obj.U_player) 
                error('U_player is empty, id: %d', obj.id);
            end
        end

        function construct_global_ineq(obj)
            
            % Construct Global Inequalities
            nc = size(obj.F_global, 1);
            
            Fh0 = [obj.F_global; zeros(obj.H * nc, obj.nx)];
            Fh = [];
            Gh = obj.G_global;
            bh = obj.alpha*obj.b_global;

            for i = 1:obj.H-1
                Fh = blkdiag(Fh, obj.F_global);
                Gh = blkdiag(Gh, obj.G_global);
                bh = [bh; obj.alpha*obj.b_global];
            end

            % Terminal Constraint
            Fh = blkdiag(Fh, obj.F_global + obj.G_global * obj.K);
            Gh = [Gh; zeros(nc, obj.H * obj.nu)];
            bh = [bh; obj.alpha*obj.b_global];

            % Final Assembly
            Fh = [zeros(nc, obj.nx * obj.H); Fh];

            Ga_x0 = Fh0 + Fh * obj.Ah; 
            Ga_u = Fh * obj.Bh + Gh;

            obj.Ga = [Ga_x0 Ga_u];
            obj.Gb = bh; 
        end

        function compute_step_sizes(obj)
            delta = 200; % Hardcoded, design parameter

            A1 = max(sum(abs(obj.Ga), 1));
            A2 = max(sum(abs(obj.Ga), 2));
            obj.TAU = 0.75 / (A1 + delta);
            obj.SIGMA = 0.75 / (A2 + 2 + delta); % Assuming weights are 1/N for each of the players. 
            obj.NU = 0.75 / (2 + delta); % Assuming weights are 1/N for each of the players. 
            fprintf('Choosing Step Sizes, TAU: %.6e , NU: %.6e, SIGMA: %.6e for player %d\n', obj.TAU, obj.NU, obj.SIGMA, obj.id);
        end

        % To get Psuedo Gradient 
        function pseudo_grad = oracle(obj, N, players)
            n_states = obj.nx + obj.H * obj.nu; 
            
            % del J / del x0 
            linear_slope_x0 = 2 * obj.state(obj.nx+1:n_states)' * obj.Bh' * obj.Qh * obj.Ah;
            quadratic_slope_x0 = 2 * obj.state(1:obj.nx)' * (obj.Q + obj.Ah' * obj.Qh * obj.Ah);

            delJ_delx0 = linear_slope_x0 + quadratic_slope_x0; 

            % del J / del U 
            linear_slope_U = 2 * obj.state(1:obj.nx)' * obj.Ah' * obj.Qh * obj.Bh + obj.fh'; 
            quadratic_slope_U = 2 * obj.state(obj.nx + 1:n_states)' * (obj.Bh' * obj.Qh * obj.Bh + obj.Rh); 
            coupled_slope_U = zeros(1, obj.H*obj.nu);
            for i = 1:N 
                if i == obj.id 
                    continue; % already considered in quadratic slope  
                end
                [player_state, ~, ~] = players{i}.communicate_state();
                coupled_slope_U = coupled_slope_U + obj.mu*(player_state(obj.nx + 1:n_states))' * obj.Rh; 
            end
            delJ_delU = linear_slope_U + quadratic_slope_U + coupled_slope_U;
            pseudo_grad = [delJ_delx0' ; delJ_delU'];
        end
    end
end