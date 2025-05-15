%% Cyclic Voltammetry Simulation
% Simulates CV responses for electron transfer reactions (E mechanism)

%% 1. System Parameters Setup

% Electrochemical Parameters
Electrode.Reaction = struct(...
    'C_red',    2,  ...    % Reduced species concentration [mmol/L]
    'C_ox',     0,  ...    % Oxidized species concentration [mmol/L]
    'k0',       1e-4, ...  % Standard rate constant [m/s]
    'E0',       0,    ...  % Formal potential [V]
    'n',        1     ...  % Electron transfer number
    );

% Electrical Parameters
Electrode.Electrical = struct(...
    'E_start',  -0.4, ...  % Initial potential [V]
    'E_switch', 0.4, ...   % Switching potential [V]
    'scan_rate', 0.02,  ...   % Scan rate [V/s]
    'num_cycles',1   ...   % Number of cycles
    );

% Physical Constants
Constants = struct(...
    'F', 96485, ...        % Faraday constant [C/mol]
    'R', 8.314,    ...     % Gas constant [J/(mol·K)]
    'T', 298,    ...       % Temperature [K]
    'D', 1e-9  ...         % Diffusion coefficient [m²/s]
    );
%% 2. Spatial and Temporal Meshing
[x,t] = CreateMesh(Electrode.Electrical, Constants);

%% 3. Potential Waveform Generation
potential = GeneratePotential(Electrode.Electrical, t);

%% 4. PDE Solution and Current Calculation
[sol, current,gradient_data] = SolvePDE(Electrode, Constants, x, t);

%% 5. Visualization
figure('Color','w','Position',[100 100 800 400])
plot(potential, current, 'LineWidth',3)
xlabel('Potential [V]'), ylabel('Current Density [A/m²]')
% title('Cyclic Voltammogram')
grid on
set(gca,'FontName','Arial','FontSize',50,'FontWeight','bold','LineWidth',2,'Box','off');
%% mesh
figure('Color','w')
plot(x, zeros(size(x)), 'bo', 'MarkerFaceColor','b')
% title('Spatial Mesh Distribution')
xlabel('Distance from Electrode (m)')
ylabel('Position')
grid on
% xlim([0 0.0001])
set(gca,'YTick',[])
set(gca,'FontName','Arial','FontSize',50,'FontWeight','bold','LineWidth',2,'Box','off');
%% Plot gradient evolution
figure;
surf(x, t, gradient_data);
xlabel('Distance (m)'); ylabel('Time (s)'); zlabel('∂[Red]/∂x');
title('Concentration Gradient Evolution');
set(gca,'FontName','Arial','FontSize',50,'FontWeight','bold','LineWidth',2,'Box','off');
%% Concentraion profile:
%% Concentration Profiles Visualization with Custom Colormaps
cRed_profiles = sol(:,:,1);
cOx_profiles = sol(:,:,2);
time_indices = [1, round(length(t)/4), round(length(t)/2), length(t)];

% Create warm (reds) and cold (blues) colormaps
num_times = length(time_indices);
warm_colors = autumn(num_times);  % Red -> Orange -> Yellow
cold_colors = winter(num_times);   % Blue -> Green

% Reverse colormaps for better temporal progression
warm_colors = flipud(warm_colors);
cold_colors = flipud(cold_colors);

figure;
hold on

% Plot concentration profiles
for idx = 1:length(time_indices)
    % Red species (warm colors)
    plot(x, cRed_profiles(time_indices(idx),:),...
        'LineWidth', 3.5,...
        'Color', warm_colors(idx,:),...
        'DisplayName', sprintf('Red t=%.1f s', t(time_indices(idx))));

    % Ox species (cold colors)
    plot(x, cOx_profiles(time_indices(idx),:),...
        'LineWidth', 3.5,...
        'Color', cold_colors(idx,:),...
        'DisplayName', sprintf('Ox t=%.1f s', t(time_indices(idx))) );
end

xlabel('Distance from Electrode (m)')
ylabel('Concentration (mmol/L)')
% title('Concentration Profiles')
legend('Location', 'northeast')
% legend off
grid on

set(gca,'FontName','Arial','FontSize',50,'FontWeight','bold','LineWidth',2,'Box','off');
% set(gca,'YScale','log')
%% Core Functions --------------------------------------------------------
function [x,t] = CreateMesh(E, C)
% Creates non-uniform spatial mesh and time vector
% Can be expressed from random walk
% L = 6*sqrt(2*abs(E.E_start-E.E_switch)/E.scan_rate*C.D);
L = 0.001;
x_near = linspace(0, L/30, 30);        % Dense near surface
x_far = linspace(1.1*L/30, L, 30);     % Sparse in bulk
x = [x_near, x_far];
% Time vector: covers complete scan cycles
t_total = 2*abs(E.E_switch-E.E_start)/E.scan_rate*E.num_cycles;
t = linspace(0, t_total, 200*E.num_cycles);
end

function potential = GeneratePotential(E, t)
% Generates triangular potential waveform
A = E.E_switch - E.E_start;
potential = sawtooth(t*pi/abs(A)*E.scan_rate, 0.5)*A/2 + (A/2 + E.E_start);
end

function [sol, current,gradient_data] = SolvePDE(E, C, x, t)
% Solves diffusion-reaction PDE system using pdepe
m = 0; % Cartesian coordinates
options = odeset('AbsTol', 1e-30);
sol = pdepe(m, @pdefun, @icfun, @bcfun, x, t, options);

% Current calculation from concentration gradient
current = zeros(size(t));
gradient_data = zeros(length(t), length(x));
for k = 1:length(t)
    % Full spatial gradient
    [~, grad] = pdeval(0, x, sol(k,:,1), x);
    gradient_data(k,:) = grad;
    % x=0 gradient
    current(k) = E.Reaction.n * C.F * C.D * grad(1);
end

%% PDE System Components
    function [c,f,s] = pdefun(~,~,u,dudx)
        % PDE: Fick's second law for both species
        c = [1; 1];         % Time derivatives
        f = C.D*[1; 1].*dudx; % Diffusion fluxes
        s = [0; 0];         % No source terms
    end

    function u0 = icfun(~)
        % Initial conditions (uniform concentrations)
        u0 = [E.Reaction.C_red; E.Reaction.C_ox];
    end

    function [pl,ql,pr,qr] = bcfun(~,ul,~,ur,t)
        % Current potential calculation
        A = E.Electrical.E_switch - E.Electrical.E_start;
        current_potential = sawtooth(t*pi/abs(A)*...
            E.Electrical.scan_rate ...
            , 0.5)*A/2 + (A/2 + E.Electrical.E_start);

        eta = current_potential - E.Reaction.E0;



        % Butler-Volmer kinetics
        alpha = 0.5;
        k0 = E.Reaction.k0;
        Fa = E.Reaction.n*C.F*alpha;
        exp_term = exp(Fa*eta/(C.R*C.T)) ;
        flux = k0*(ul(1)*exp_term - ul(2)/exp_term);



        % Boundary conditions
        pl = [flux; flux];  % Proper flux terms for both species
        ql = [-1; 1];      % Flux directions
        pr = [0; 0];       % No source terms
        qr = [1; 1];       % Zero flux conditions
    end

end

