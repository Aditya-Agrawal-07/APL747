close all;
clear all;

% Model parameters
eps = 1;
sigma = 1;
mass = 1;
L = 6.8;
kb = 1.38e-23;

% Reading positions from liquid256.txt
r = load("liquid256.txt");

N = size(r,1); dim = size(r,2);

% Initialization of velocities
v = -1 + 2*rand(N,dim);

% Net momentum
p_net = mass*sum(v);

% Adjusted velocities
v = v - p_net/(mass*N);

% Updated net momentum 
p_net = mass*sum(v);

% Timestepping parameters
timesteps = 1e5;
final_time = 200;
dt = final_time/timesteps;
time = (0:timesteps)*dt; % Time vector

% Arrays to store simulation variables
KE = zeros(timesteps+1,1);
PE = zeros(timesteps+1,1);
TE = zeros(timesteps+1,1);
T_inst = zeros(timesteps+1,1);
p = zeros(timesteps+1,dim);
f = zeros(timesteps+1,N,dim);

% Initial values of parameters
KE(1) = Kinst(v);
y = FVcalculator(r);
PE(1) = y(size(y,1),1);
TE(1) = KE(1) + PE(1);
T_inst(1) = 2*KE(1)/(3*N*kb);
p(1,:) = p_net(:);

% Initial values of force, position, velocity
f_t = y(1:N,:);
f(1,:,:) = f_t(:,:);
r_t = r;
v_t = v;

% Timestepping using Velocity Verlet Scheme
for t=1:timesteps

    disp(t);

    r_t = r_t + dt*v_t + (dt^2/(2*mass))*f_t;
    % Replacing position of atoms that move out of the box with their
    % images
    for i=1:N
        for k=1:dim
            if r_t(i,k) > L
                r_t(i,k) = r_t(i,k) - L;
            elseif r_t(i,k) < 0
                r_t(i,k) = r_t(i,k) + L;
            end
        end
    end
    
    f0 = f_t;
    y = FVcalculator(r_t);
    f_t = y(1:N,:);
    v_t = v_t + (dt/(2*mass))*(f_t+f0);

    KE(t+1) = Kinst(v_t);
    PE(t+1) = y(N+1,1);
    TE(t+1) = KE(t+1) + PE(t+1);
    p(t+1,:) = mass*sum(v_t);
    T_inst(t+1) = 2*KE(t+1)/(3*N*kb);
    f(t+1,:,:) = f_t(:,:);

end

disp(mean(T_inst));
disp(mean(KE));

% Plotting KE, PE, and TE
figure;
plot(time, KE, 'r', 'LineWidth', 1.5); 
hold on;
plot(time, PE, 'b', 'LineWidth', 1.5); 
plot(time, TE, 'g', 'LineWidth', 1.5); 
hold off;
xlabel('Time (LJ units)');
ylabel('Energy (\epsilon)');
title('Kinetic, Potential, and Total Energy vs Time');
legend('Kinetic Energy', 'Potential Energy', 'Total Energy');
grid on;

% Plotting KE
figure;
plot(time, KE, 'b', 'LineWidth', 1.5);
xlabel('Time (LJ units)');
ylabel('Energy (\epsilon)');
title('Kinetic Energy Vs Time');
grid on;

% Plotting PE
figure;
plot(time, PE, 'b', 'LineWidth', 1.5);
xlabel('Time (LJ units)');
ylabel('Energy (\epsilon)');
title('Potential Energy Vs Time');
grid on;

% Plotting TE
figure;
plot(time, TE, 'b', 'LineWidth', 1.5);
xlabel('Time (LJ units)');
ylabel('Energy (\epsilon)');
title('Total Energy Vs Time');
grid on;

% Plotting Instantaneous Temperature (T_inst)
figure;
plot(time, T_inst, 'b', 'LineWidth', 1.5); 
xlabel('Time (LJ units)');
ylabel('Temperature (\epsilon/kB)');
title('Instantaneous Temperature vs Time');
grid on;

% Plotting momentum components (Px, Py, Pz)
figure;
plot(time, p(:,1), 'b', 'LineWidth', 1.5);
xlabel('Time (LJ units)');
ylabel('P_x (LJ units)');
title('Momentum in X direction vs Time');
grid on; 

figure;
plot(time, p(:,2), 'b', 'LineWidth', 1.5); 
xlabel('Time (LJ units)');
ylabel('P_y (LJ units)');
title('Momentum in Y direction vs Time');
grid on; 

figure;
plot(time, p(:,3), 'b', 'LineWidth', 1.5);
xlabel('Time (LJ units)');
ylabel('P_z (LJ units)');
title('Momentum in Z direction vs Time');
grid on;

% Function to evaluate force vector and total Potential
function y = FVcalculator(r)

    % Model parameters
    eps = 1;
    sigma = 1;
    L = 6.8;
    N = size(r,1); dim = size(r,2);

    % Calculation of interatomic potential energy and force
    % and implementation of cutoff radius scheme
    V = 0;
    F = zeros(N,dim);
    r_cut = 2.5;
    f_r_cut = (24/r_cut^7)*((2/r_cut^6)-1);
    V_r_cut = (4/r_cut^6)*((1/r_cut^6)-1);
    
    for i=1:N-1
        for j=i+1:N
    
            rij = r(i,:) - r(j,:);
            % Implementation of periodic B.C.
            for k=1:dim
                if rij(k) > L/2
                    rij(k) = rij(k) - L;
                elseif rij(k) < -L/2
                    rij(k) = rij(k) + L;
                end
            end
            r_ij = norm(rij);
            
            % Implementation of cutoff scheme
            if r_ij <= r_cut
                V_int = (4/r_ij^6)*((1/r_ij^6) - 1) - V_r_cut + (r_ij-r_cut)*f_r_cut;
                f_int = (24/r_ij^7)*((2/r_ij^6) - 1) - f_r_cut;
        
                F(i,:) = F(i,:) + f_int*rij(1,:)/r_ij;
                F(j,:) = F(j,:) - f_int*rij(1,:)/r_ij;
    
                V = V + V_int;
            end
        end
    end 
    V_ = [V,0,0];
    y = [F; V_];
end

% function to calculate Kinetic Energy
function y = Kinst(v)

    % Instantaneous Kinetic Energy
    K_inst = 0;
    n = size(v,1);
    mass = 1;
    for p=1:n
        K_inst = K_inst + 0.5*mass*norm(v(p,:))^2;
    end

    y = K_inst;
end


