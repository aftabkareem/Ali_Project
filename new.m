% Parameters
L = 1;              % Length of the domain
Nx = 200;           % Number of grid points
dx = L / Nx;        % Grid spacing
x = linspace(0, L, Nx); % Spatial grid
dt = 0.01;          % Time step
T = 0.3;            % Total simulation time
Nt = ceil(T / dt);  % Number of time steps
u_values = [0.25, 0.4, 0.5]; % Different velocity values

% Initialize solution matrix
u_solution = zeros(Nx, Nt, length(u_values));

% Loop over different velocity values
for u_index = 1:length(u_values)
    u = u_values(u_index);
    
    % Stabilized discretization coefficient
    c = u * dt / dx;
    
    % Check CFL condition
    if c > 1
        error('CFL condition violated: 0 <= c <= 1');
    end
    
    % Initial condition (Gaussian wave with width 0.05)
    sigma = 0.05; % Width of the Gaussian wave
    u_solution(:, 1, u_index) = exp(-((x-0.5*L).^2) / (2*sigma^2));
    
    % Time integration using implicit upwind scheme
    A = sparse(Nx, Nx);
    for i = 2:Nx
        A(i, i) = 1 + c;
        A(i, i-1) = -c;
    end
    A(1, 1) = 1 + c;
    A(1, Nx) = -c;
    
    for n = 1:Nt-1
        % Solve the implicit equation
        u_solution(:, n+1, u_index) = A \ u_solution(:, n, u_index);
    end
    
    % Plot solution
    figure;
    for n = 1:10:Nt
        plot(x, u_solution(:, n, u_index));
        hold on;
    end
    hold off;
    title(['Solution for CFL = ', num2str(c)]);
    xlabel('x');
    ylabel('\Phi');
    grid on;
    legend('t = 0', 't = 0.01', 't = 0.02');
end
