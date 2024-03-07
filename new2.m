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
    
    % Time integration
    for n = 1:Nt-1
        % Apply stabilized discretization
        for i = 2:Nx
            u_solution(i, n+1, u_index) = u_solution(i, n, u_index) - (u * dt / dx) * (u_solution(i, n, u_index) - u_solution(i-1, n, u_index));
        end
        
        % Boundary condition (periodic)
        u_solution(1, n+1, u_index) = u_solution(Nx, n+1, u_index);
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
