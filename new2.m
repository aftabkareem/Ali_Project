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



% new-------------new

% Parameters
L = 1;                  % Length of the domain
Nx_values = [100, 200, 400]; % Different number of grid points for testing
T = 3;                  % Total simulation time
u_values = [0.25, 0.4, 0.5]; % Different velocity values
times_to_plot = [0, 0.15, 0.45, 0.65, 0.75]; % Times at which we want to plot the solution

% Loop over different number of grid points
for Nx_index = 1:length(Nx_values)
    Nx = Nx_values(Nx_index);
    dx = L / Nx;        % Grid spacing
    x = linspace(0, L, Nx); % Spatial grid
    dt = 0.01;          % Time step
    Nt = ceil(T / dt);  % Number of time steps

    % Initialize solution matrix
    u_solution = zeros(Nx, Nt, length(u_values));
    
    % Measure the speed of calculation
    tic;
    
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
            % Apply upwind discretization with periodic boundary condition
            u_solution(1, n+1, u_index) = u_solution(1, n, u_index) - c * (u_solution(1, n, u_index) - u_solution(Nx, n, u_index));
            for i = 2:Nx
                if u >= 0
                    % Upwind scheme if velocity is positive
                    u_solution(i, n+1, u_index) = u_solution(i, n, u_index) - c * (u_solution(i, n, u_index) - u_solution(i-1, n, u_index));
                else
                    % Upwind scheme if velocity is negative
                    u_solution(i, n+1, u_index) = u_solution(i, n, u_index) - c * (u_solution(i+1, n, u_index) - u_solution(i, n, u_index));
                end
            end
        end
        
        % Determine time for wave to leave the domain
        wave_exit_time = find_wave_exit_time(u_solution(:,:,u_index), dt);
        
        % Capture locations of the wave at various time instances
        wave_locations = capture_wave_locations(u_solution(:,:,u_index), times_to_plot, dt);
        
        % Plot solution for specific times
        figure;
        legends = {};
        for t_index = 1:length(times_to_plot)
            t_value = times_to_plot(t_index);
            n = round(t_value / dt) + 1; % Calculate the time step (adding 1 because MATLAB indexing starts at 1)
            plot(x, u_solution(:, n, u_index));
            hold on;
            legends{end+1} = ['Time = ' num2str(t_value, '%.2f') ' s, Location = ' num2str(wave_locations(t_index))];
        end
        hold off;
        
        title(['Solution for u = ', num2str(u), ', Nx = ', num2str(Nx), ', CFL = ', num2str(c)]);
        xlabel('x');
        ylabel('\Phi');
        grid on;
        legend(legends, 'Location', 'best');
        fprintf('Wave exit time for u = %.2f, Nx = %d: %.2f s\n', u, Nx, wave_exit_time);
    end
    
    % End measurement of speed
    computation_time = toc;
    fprintf('Computation time for Nx = %d: %.2f s\n', Nx, computation_time);
end

% Function to determine time for wave to leave the domain
function exit_time = find_wave_exit_time(u_solution, dt)
    exit_time = NaN;
    threshold = max(u_solution(:,1)) * 0.01; % 1% of the initial peak value
    for n = 1:size(u_solution, 2)
        if max(u_solution(:,n)) < threshold
            exit_time = (n-1) * dt;
            return;
        end
    end
end

% Function to capture locations of the wave's peak at various time instances
function locations = capture_wave_locations(u_solution, times_to_plot, dt)
    locations = zeros(1, length(times_to_plot));
    for t_index = 1:length(times_to_plot)
        n = round(times_to_plot(t_index) / dt) + 1;
        [~, loc] = max(u_solution(:, n));
        locations(t_index) = loc;
    end
end

