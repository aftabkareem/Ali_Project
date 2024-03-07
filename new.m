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




--new---

% Parameters
L = 1;              % Length of the domain
Nx_values = [50, 100, 200]; % Different numbers of grid points for comparison
dx_values = L ./ Nx_values; % Grid spacing for each Nx value
dt = 0.01;          % Time step
T = 0.45;            % Total simulation time
Nt = ceil(T / dt);  % Number of time steps
u_value = 0.25;      % Wave velocity (chosen to satisfy CFL)

% Initialize a matrix to store time for the wave to leave the domain for different Nx
time_to_leave = zeros(length(Nx_values), 1);

% Loop over different numbers of nodes
for idx = 1:length(Nx_values)
    Nx = Nx_values(idx);
    dx = dx_values(idx);
    x = linspace(0, L, Nx); % Spatial grid
    
    % Stabilized discretization coefficient
    c = u_value * dt / dx;
    
    % Check CFL condition
    if c > 1
        error('CFL condition violated for u = %f: c (%f) must be <= 1', u_value, c);
    end
    
    % Initial condition (Gaussian wave with width 0.05)
    sigma = 0.05; % Width of the Gaussian wave
    u_initial = exp(-((x-0.5*L).^2) / (2*sigma^2));
    
    % Time integration using implicit upwind scheme
    A = spdiags([c*ones(Nx, 1), (1-c)*ones(Nx, 1)], [-1, 0], Nx, Nx);
    A(1, Nx) = c; % Periodic boundary condition
    
    % Initialize solution array
    u_solution = zeros(Nx, Nt);
    u_solution(:, 1) = u_initial;
    
    tic; % Start timing for speed of calculation
    for n = 1:Nt-1
        % Solve the implicit equation
        u_solution(:, n+1) = A \ u_solution(:, n);
    end
    computation_time = toc; % Time for calculation
    
    % Find the time when the peak of the wave first reaches the end of the domain
    [~, max_indices] = max(u_solution, [], 1); % Find peak indices for each time step
    wave_exits = find(x(max_indices) >= L, 1); % Find the first time the peak exits the domain
    
    if ~isempty(wave_exits)
        time_to_leave(idx) = wave_exits * dt;
    else
        time_to_leave(idx) = NaN; % If wave never exits during the simulation
    end
    
    % Display speed of calculation
    disp(['Computation time for Nx = ', num2str(Nx), ': ', num2str(computation_time), ' seconds']);
    
    % Capture the locations of the wave at t=0, 0.15, and 0.45
    times_to_capture = [0, 0.15, 0.45];
    figure;
    for capture_time = times_to_capture
        time_index = ceil(capture_time / dt) + 1;
        plot(x, u_solution(:, time_index));
        hold on;
    end
    title(['Wave propagation for Nx = ', num2str(Nx)]);
    xlabel('Domain (x)');
    ylabel('Wave amplitude (\Phi)');
    legend(arrayfun(@(t) ['t = ', num2str(t)], times_to_capture, 'UniformOutput', false));
    hold off;
end

% Display time for the wave to leave the domain for different Nx
disp('Time for the wave to leave the domain for different Nx values:');
disp(time_to_leave);
