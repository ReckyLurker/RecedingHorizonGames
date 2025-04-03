% Simulates H steps of the game. 
clear all; 
clc; 

%%%%%%%%% Add Mosek Path [Possibly the Below] %%%%%%%%%%%%
% addpath("C:\Program Files\Mosek\11.0\tools\platform\win64x86\bin")
% addpath("C:\Program Files\Mosek\11.0\toolbox\r2019b")

%{
   If Mosek is not available, use Sedumi-1.3. It comes with MPT3 Toolbox. Change Line 325 of Agent.m file from 'mosek' to 'sedumi' in the project_point routine. 
   Run this code from project root. 
%}

% Import MPT3 
addpath('libraries/mpt3/tbxmanager/')
addpath('receding-horizon/src/utils/')
addpath('receding-horizon/src/systems/')
startup % Run startup file to add all required libraries to the path.
import mpt.*

%%%%%%%%%%%%%%%%%%%%%%%%%%% GAME SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Number of Players and Horizon Length %%%%
N = 4; % Number of Active Players. They consume and also store in batteries.
H = 10; % Horizon Length for the Game-Theoretic MPC
MAX_ITER = 7500; % Maximum number of iterations to find vGNE 
MIN_ITER = 100; % Minimum number of iterations to find vGNE 
T = 50; % Total Time-Steps to simulate
MONTE_CARLO_SIMULATIONS = 200; 

%%%% Player Dynamics Setup %%%%
nx = 2; % Number of States 
nu = 2; % Number of Control Inputs 

%%% Shared Constraints %%%
L_max = 185; % Maximum Grid Capcity 
L_min = 0; % Minimum Grid Capacity 

%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%
players = {};
q_ref = {};
e_ref_center = {};
x = {};
for i = 1:N 
   if i == 1 
      A = [1 0; 0 0.9956];
      B = [1 0; 0 0.9];
      Q = 10*[0.0655 0; 0 0.19];
      R = 0.0413*ones(nu,nu);
      f = 4*0.0413*ones(nu,1);
      mu = 1; 
      q_ref{i} = 0.6*55;
      x{i} = [0; -q_ref{i}*0.9];
      e_ref_center{i} = 0.5*(25 + 35);
      X = state_polyhedron([-7; 7], [0; 55], q_ref{i});
      U = control_polyhedron([0; 37], [-5.5; 5.5], [0; 41], B(2,2), e_ref_center{i});
      W = disturbance_polyhedron([25; 35], A(2,2), q_ref{i});
   elseif i == 2
      A = [1 0; 0 0.98];
      B = [1 0; 0 0.8882];
      Q = 10*[0.0193 0; 0 0.32];
      R = 0.0236*ones(nu,nu);
      f = 4*0.0236*ones(nu,1);
      mu = 1;
      q_ref{i} = 0.5*65; 
      x{i} = [-10; -0.9*q_ref{i}];
      e_ref_center{i} = 0.5*(30 + 45);
      X = state_polyhedron([-10; 10], [0; 65], q_ref{i});
      U = control_polyhedron([0; 47], [-6.5; 6.5], [0; 51], B(2,2), e_ref_center{i});
      W = disturbance_polyhedron([30; 45], A(2,2), q_ref{i});
   elseif i == 3
      A = [1 0; 0 0.96];
      B = [1 0; 0 0.9084];
      Q = 10*[0.0445 0; 0 0.6];
      R = 0.0342*ones(nu,nu);
      f = 4*0.0342*ones(nu,1);
      mu = 1;
      q_ref{i} = 0.6*35;
      x{i} = [5; -q_ref{i}*0.9];
      e_ref_center{i} = 0.5*(15+32);
      X = state_polyhedron([-11; 11], [0; 35], q_ref{i});
      U = control_polyhedron([0; 34], [-3.5; 3.5], [0; 36], B(2,2), e_ref_center{i});
      W = disturbance_polyhedron([15; 32], A(2,2), q_ref{i});
   else 
      A = [1 0; 0 0.94];
      B = [1 0; 0 0.9860];
      Q = 10*[0.0232 0; 0 0.4];
      R = 0.02765*ones(nu,nu);
      f = 4*0.02765*ones(nu,1);
      mu = 1;
      q_ref{i} = 0.6*45; 
      x{i} = [-7; -q_ref{i}*0.9];
      e_ref_center{i} = 0.5*(22+42);
      X = state_polyhedron([-13; 13], [0; 45], q_ref{i});
      U = control_polyhedron([0; 45], [-6.5; 6.5], [0; 45], B(2,2), e_ref_center{i});
      W = disturbance_polyhedron([22; 42], A(2,2), q_ref{i});
   end 
   % Players share similar constraints. 
   F_global = [0 0; 0 0];
   G_global = [-1 -1; 1 1];
   b_global = [-(L_min - M * p_max - 123); (L_max - M * p_max - 123)];
   players{i} = Agent(i, H, A, B, Q, R, f, mu, X, U, W, x{i}, F_global, G_global, b_global, 1/N);
end

% Tighten Shared Constraints 
for i = 1:N 
   for j = 1:N 
      delta = players{j}.communicate_shared_constraint_tightening();
      players{i}.update_shared_constraint(delta);
   end
end

% Prepare Steady State Game (Cost-Function, Constraints, etc.)
for i = 1:N 
   players{i}.compute_steady_state_parameters();
end

% Compute the steady state vGNE 
iter = 0;
tic; 
while iter < MAX_ITER 
   for j = 1:N 
      players{j}.update_ss_state(N, players); 
   end

   for j = 1:N 
      players{j}.update_ss_lambda(N, players);
   end

   state_diff = 0;
   lamda_diff = 0;
   s = 0;
   for j = 1:N 
      state_diff = state_diff + (1/N) * norm(players{j}.ss_tstate - players{j}.ss_state) / (norm(players{j}.ss_state) + 1e-4);
      for k = j+1:N 
         [~, new_lamdaj, ~] = players{j}.communicate_ss_update();
         [~, new_lamdak, ~] = players{k}.communicate_ss_update(); 
         lamda_diff = lamda_diff + norm(players{k}.ss_tlambda - players{j}.ss_tlambda) / (norm(players{j}.ss_tlambda) + 1e-4);
         s = s + 1;
      end
   end
   lamda_diff = lamda_diff / s;
   for j = 1:N 
      players{j}.finish_ss_update();
   end
   iter = iter + 1;
   % fprintf('Curr Iter: %d iter with consensus: %.6e, state_diff: %.6e\n', iter, lamda_diff, state_diff);
   if 0.5*(state_diff + lamda_diff) < 1e-4 && iter > MIN_ITER 
      break; 
   end
end
time_taken = toc; 
fprintf('Steady State vGNE solution found after %d iter with consensus: %.6e, state_diff: %.6e\n', iter, lamda_diff, state_diff);
fprintf('Steady State vGNE Loop Exec Time: %.6f\n', time_taken);

for i = 1:N 
   fprintf('Steady State for Player %d:\n', i);
   disp(players{i}.ss_state);
   fprintf('\n\n\n');
end

% Shift the system 
for i = 1:N 
   players{i}.shift_system(); % Shift Local Constraints 
   for j = 1:N 
      delta = players{j}.F_global * players{j}.ss_state(1:players{j}.nx) + players{j}.G_global * players{j}.ss_state(players{j}.nx+1:players{j}.nx+players{j}.nu);
      players{i}.update_shared_constraint(delta); % Shift Shared Constraints 
   end
   players{i}.construct_terminal_set(); % Compute Terminal Constraint
   players{i}.prepare_sim(); % Prepare the Local Polyhedron (without dissipativity), (Global Inequalities), (Step Sizes for Optimization)
end

% Store Simulation Data
zeta_over_time = {};
q_over_time = {};
e_over_time = {};
e_ref_over_time = {};
s_over_time = {};
net_load_over_time = [];

x0 = {};
for i = 1:N 
   x0{i} = x{i}; 
end

% Store Monte-Carlo data 
zeta_over_time_MC = {};
q_over_time_MC = {};
e_over_time_MC = {};
e_ref_over_time_MC = {};
s_over_time_MC = {};
net_load_over_time_MC = [];

for i = 1:N 
   zeta_over_time_MC{i} = [];
   q_over_time_MC{i} = [];
   e_over_time_MC{i} = [];
   e_ref_over_time_MC{i} = [];
   s_over_time_MC{i} = [];
end

for progress = 1:MONTE_CARLO_SIMULATIONS 
   tic; 
   for i = 1:N 
      zeta_over_time{i} = [x0{i}(1)];
      q_over_time{i} = [0];
      e_over_time{i} = [];
      e_ref_over_time{i} = [];
      s_over_time{i} = [];
      net_load_over_time = zeros(T, 1);
      players{i}.update_polyhedron(x0{i});
      x{i} = x0{i};
   end
   
   for t = 1:T 
      iter = 0; 
      % tic; 
      while iter < MAX_ITER 
         for i = 1:N 
            players{i}.update_state(N, players, t);
         end 
         for i = 1:N 
            players{i}.update_lambda(N, players);
         end
         state_diff = 0;
         lamda_diff = 0;
         s = 0;
         for i = 1:N 
            state_diff = state_diff + (1/N) * norm(players{i}.tstate - players{i}.state) / (norm(players{i}.state) + 1e-6);
            for j = i+1:N 
               [~, new_lamdaj, ~] = players{i}.communicate_update();
               [~, new_lamdak, ~] = players{j}.communicate_update(); 
               lamda_diff = lamda_diff + norm(players{j}.tlambda - players{i}.tlambda) / (norm(players{i}.tlambda) + 1e-6);
               s = s + 1;
            end
         end
         lamda_diff = lamda_diff / s;
         for i = 1:N 
            players{i}.finish_update();
         end
         iter = iter + 1;
         % fprintf('Curr Iter: %d iter with consensus: %.6e, state_diff: %.6e\n', iter, lamda_diff, state_diff);
         if 0.5*(state_diff + lamda_diff) < 1e-4 && iter > MIN_ITER 
            break; 
         end
      end
      % time_taken = toc; 
      % fprintf('vGNE Loop Exec Time: %.6f\n', time_taken);
      fprintf('vGNE solution found after %d iter with consensus: %.6e, state_diff: %.6e\n', iter, lamda_diff, state_diff);
   
      for i = 1:N 
         [x{i}, u, e_ref_t] = players{i}.next(x{i});
         [J_player, J_stage] = players{i}.predict_cost(N, players);
         players{i}.update_polyhedron(x{i});
         fprintf('Player %d Nominal Stage Cost: %.6f\n', i, J_stage);
         fprintf('Player %d Cost: %.6e\n', i, J_player);
   
         zeta_over_time{i} = [zeta_over_time{i}; x{i}(1)];
         q_over_time{i} = [q_over_time{i}; x{i}(2) + q_ref{i}];
         e_over_time{i} = [e_over_time{i}; u(1) + e_ref_center{i}];
         s_over_time{i} = [s_over_time{i}; u(2)];
         e_ref_over_time{i} = [e_ref_over_time{i}; -e_ref_t(1) + e_ref_center{i}];
         net_load_over_time(t) = net_load_over_time(t) + u(1) + e_ref_center{i} + u(2); 
         players{i}.update_SDP_constraint(); % Update Dissipativity constraint 
      end
      fprintf('Current Iter: %d\n', t);
   end
   fprintf('Current MONTE_CARLO_ITERATION: %d\n', progress);
   % figure_gen; 
   for i = 1:N 
      zeta_over_time_MC{i} = [zeta_over_time_MC{i} zeta_over_time{i}];
      q_over_time_MC{i} = [q_over_time_MC{i} q_over_time{i}];
      e_over_time_MC{i} = [e_over_time_MC{i} e_over_time{i}];
      s_over_time_MC{i} = [s_over_time_MC{i} s_over_time{i}];
      e_ref_over_time_MC{i} = [e_ref_over_time_MC{i} e_ref_over_time{i}];
   end
   net_load_over_time_MC = [net_load_over_time_MC net_load_over_time];
   time_taken = toc; 
   fprintf('\n\nCurrent MC Iteration took: %.6f seconds\n\n', time_taken);
end

% Store to CSV 
for i = 1:N 
   csvwrite(sprintf('simulation_data/MC_zeta_%d.csv', i), zeta_over_time_MC{i});
   csvwrite(sprintf('simulation_data/MC_q_%d.csv', i), q_over_time_MC{i});
   csvwrite(sprintf('simulation_data/MC_e_%d.csv', i), e_over_time_MC{i});
   csvwrite(sprintf('simulation_data/MC_s_%d.csv', i), s_over_time_MC{i});
   csvwrite(sprintf('simulation_data/MC_e_ref_%d.csv', i), e_ref_over_time_MC{i});
end
csvwrite('simulation_data/MC_netload.csv', net_load_over_time_MC);

% Convinience Functions 
function X = state_polyhedron(zeta_constraint, q_constraint, q_ref)
   x_min = [zeta_constraint(1); q_constraint(1) - q_ref];
   x_max = [zeta_constraint(2); q_constraint(2) - q_ref];
   X = Polyhedron('lb', x_min, 'ub', x_max); % State Polyhedron
end

function U = control_polyhedron(e_constraint, s_constraint, load_constraint, beta, e_ref_center)
   U_A = [-1 0; 0 -beta; 1 0; 0 beta; -1 -1; 1 1];
   U_b = [
      -e_constraint(1) + e_ref_center;
      -s_constraint(1);
      e_constraint(2) - e_ref_center; 
      s_constraint(2);
      -load_constraint(1) + e_ref_center;
      load_constraint(2) - e_ref_center;
   ];
   U = Polyhedron(U_A, U_b); 
end

function W = disturbance_polyhedron(e_ref, alpha, q_ref)
   e_ref_center = (e_ref(1) + e_ref(2)) / 2; 
   w_min = [-(e_ref(2) - e_ref_center); (alpha - 1) * q_ref];
   w_max = [-(e_ref(1) - e_ref_center); (alpha - 1) * q_ref];
   W = Polyhedron('lb', w_min, 'ub', w_max);
end