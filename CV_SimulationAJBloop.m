close all;
clear;
clc;

%% Import files:   THIS VERSION OF THE CODE IS TO PLACE IN THE FILE WITH THE DATA TO ANALYSE; IT DEFAULTS TO ITS CURRENT DIRECTORY

cd    % <<<< So that the Matlab script defaults to the directory that it is put in, no other directory is specified!
directory = cd     ;

txt_data_kET = dir ("*k_ET_0_CVsim*.xlsx");
[folderPath, folderName, ~] = fileparts(directory);
filename_kET = txt_data_kET.name;

% Read the k_ET_0 m/s data from the first sheet
data_kET = readtable(filename_kET, 'Sheet', 1);

% Display the first few rows of the table
disp(data_kET);

% Extract timepoints and specific columns (example)
timepoints = data_kET{:,2}; % Assuming first column is Timepoint
kET_values = table2array(data_kET(:, 3:end)); % Extracting the rest of the data 

% Create a logical matrix: 1 for data present, 0 for missing/NaN
logicalMatrix = ~isnan(kET_values);

% Display the logical matrix
disp(logicalMatrix);
disp(sum(sum(logicalMatrix)));

% Read the D in m2/s data
txt_data_D = dir ("*D_CVsim*.xlsx");
[folderPath, folderName, ~] = fileparts(directory);
filename_D = txt_data_D.name;

% Read the data from the first sheet
data_D = readtable(filename_D, 'Sheet', 1);

% Display the first few rows of the table
disp(data_D);

% Extract timepoints and specific columns (example)
D_values = table2array(data_D(:, 3:end)); % Extracting the rest of the data

% Define the number of rows based on your data (same as logical matrix size)
numRows = size(logicalMatrix, 1); % Same number of rows as data

% Define the concentration matrix
concentrationMatrix = [1 * ones(numRows, 7), 0.2 * ones(numRows, 5)];

% Display the matrix
disp(concentrationMatrix);

% Extract all column headers (variable names) from the table
headers_kET = data_kET.Properties.VariableNames(3:end);

%% Cyclic Voltammetry Simulation
% Simulates CV responses for electron transfer reactions (E mechanism)


%% 1. System DEFAULT Parameters Setup

% Electrochemical Parameters
Electrode.Reaction = struct(...
    'C_red',    1,  ...    % Reduced species concentration [mmol/L]
    'C_ox',     0,  ...    % Oxidized species concentration [mmol/L]
    'k0',       1e-1 , ...  % Standard rate constant [m/s]
    'E0',       0.195,    ...  % Formal potential [V]
    'n',        1     ...  % Electron transfer number
    );

% Electrical Parameters
Electrode.Electrical = struct(...
    'E_start',  -0.3, ...  % Initial potential [V]
    'E_switch', 0.7, ...   % Switching potential [V]
    'scan_rate', 0.05,  ...   % Scan rate [V/s]
    'num_cycles',5   ...   % Number of cycles
    );

% Physical Constants
Electrode.Constants = struct(...
    'F', 96485, ...        % Faraday constant [C/mol]
    'R', 8.314,    ...     % Gas constant [J/(mol·K)]
    'T', 298,    ...       % Temperature [K]
    'D', 6.76876E-10  ...         % Diffusion coefficient [m²/s]           % 1.0 = An unreasonable default value (CV A-zone) --> Will be obvious if something wrong is plotted
       );                                                                  % Reasonable default: 6.76876E-10 m2/s from EIS or 1.73102E-09 m2/s from CV Randles-Sevcik

%% Now make a structure of matrices with values extracted from EIS fitting (Circ1, Randles-CPE-W)

% Get matrix dimensions
[numRows, numCols] = size(kET_values); 

% Initialize a struct array to hold the parameter sets for each condition
ElectrodeArray = repmat(Electrode, numRows, numCols);

% Loop through the matrix to assign values dynamically
for i = 1:numRows
    for j = 1:numCols
        if logicalMatrix(i, j) % Check if data exists at this point
            % Create a new struct entry
            ElectrodeArray(i, j).Reaction = struct( ...
                'C_red', concentrationMatrix(i, j), ...  % Set C_red from concentration matrix
                'C_ox', concentrationMatrix(i, j)-concentrationMatrix(i, j), ...                           % Default value
                'k0', kET_values(i, j), ...            % Set k0 from data matrix
                'E0', 0.195, ...                         % Keep other reaction parameters constant
                'n', 1 ...
            );

            % Assign diffusion coefficient D:                                      <<<<<<<<<<<<<<<<<<<<< Enable if you don't want to use standard D!
             ElectrodeArray(i, j).Constants.D = D_values(i, j);                  
        end
    end
end

% Display example parameters for the first timepoint and condition
disp(ElectrodeArray(1,1).Reaction);

%%  Plot

% Loop through each row (timepoints) and column (electrode conditions)
for j = 1:size(logicalMatrix, 2)        %column (expt)
    
    figure('Color','w','Position',[100 100 800 400])

    for i = 1:size(logicalMatrix, 1)    %row (scan/timepoint)
        
        % Check if the data exists (logicalMatrix is true)
        if logicalMatrix(i, j) == 1  
            
            % Extract parameters for the current valid entry
            Electrode = ElectrodeArray(i, j);  
            
            % 2. Spatial and Temporal Meshing
            [x, t] = CreateMesh(Electrode.Electrical, Electrode.Constants);
            
            % 3. Potential Waveform Generation
            potential = GeneratePotential(Electrode.Electrical, t);
            
            % 4. PDE Solution and Current Calculation
            [~, current, ~] = SolvePDE(Electrode, Electrode.Constants, x, t);
            
            % 5. Visualization
            
            timepoint = num2str(timepoints(i));
            fprintf('the current i and j is %d,%d\n',i,j)
            fprintf([timepoint '\n'])
            plot(potential(800:1000), current(800:1000), 'LineWidth', 3, DisplayName=num2str(timepoints(i)))
            xlabel('Potential [V]'), ylabel('Current Density [A/m²]')
            title(sprintf('Electrode Condition at (Column: %d, Row: %d)', j, i))
            grid on
            set(gca, 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold', 'LineWidth', 2, 'Box', 'off');

            legend show

            hold on;

            fprintf('%d_%d_done\n', j, i);

            % Export Cycle 2:
            exportMatrix2(:,1) = potential(200:400);
            exportMatrix2(:,2) = current(200:400);

            
            % Export Cycle 5:
            exportMatrix5(:,1) = potential(800:1000);
            exportMatrix5(:,2) = current(800:1000);

            % Create the filename using num2str
            timepoint = num2str(timepoints(i),'%.2f');
            filename = sprintf('%s_%sh_%s.txt', headers_kET{j}, timepoint, 'BVsim_cyc2')

            % Save exportMatrix to the file
            writematrix(exportMatrix2, filename, 'Delimiter', 'tab')

            % Create the filename using num2str
            filename = sprintf('%s_%sh_%s.txt', headers_kET{j}, timepoint, 'BVsim_cyc5')

            % Save exportMatrix to the file
            writematrix(exportMatrix5, filename, 'Delimiter', 'tab')

        end
    end

    hold off;

end




% 
% % 2. Spatial and Temporal Meshing
% [x,t] = CreateMesh(Electrode.Electrical, Electrode.Constants);
% 
% % 3. Potential Waveform Generation
% potential = GeneratePotential(Electrode.Electrical, t);
% 
% % 4. PDE Solution and Current Calculation
% [sol, current,gradient_data] = SolvePDE(Electrode, Electrode.Constants, x, t);
% 
% % 5. Visualization
% figure('Color','w','Position',[100 100 800 400])
% plot(potential, current, 'LineWidth',3)
% xlabel('Potential [V]'), ylabel('Current Density [A/m²]')
% % title('Cyclic Voltammogram')
% grid on
% set(gca,'FontName','Arial','FontSize',50,'FontWeight','bold','LineWidth',2,'Box','off');
% 
% 


%% mesh
% figure('Color','w')
% plot(x, zeros(size(x)), 'bo', 'MarkerFaceColor','b')
% % title('Spatial Mesh Distribution')
% xlabel('Distance from Electrode (m)')
% ylabel('Position')
% grid on
% % xlim([0 0.0001])
% set(gca,'YTick',[])
% set(gca,'FontName','Arial','FontSize',50,'FontWeight','bold','LineWidth',2,'Box','off');
% %% Plot gradient evolution
% figure;
% surf(x, t, gradient_data);
% xlabel('Distance (m)'); ylabel('Time (s)'); zlabel('∂[Red]/∂x');
% title('Concentration Gradient Evolution');
% set(gca,'FontName','Arial','FontSize',50,'FontWeight','bold','LineWidth',2,'Box','off');
% %% Concentraion profile:
% %% Concentration Profiles Visualization with Custom Colormaps
% cRed_profiles = sol(:,:,1);
% cOx_profiles = sol(:,:,2);
% time_indices = [1, round(length(t)/4), round(length(t)/2), length(t)];
% 
% % Create warm (reds) and cold (blues) colormaps
% num_times = length(time_indices);
% warm_colors = autumn(num_times);  % Red -> Orange -> Yellow
% cold_colors = winter(num_times);   % Blue -> Green
% 
% % Reverse colormaps for better temporal progression
% warm_colors = flipud(warm_colors);
% cold_colors = flipud(cold_colors);
% 
% figure;
% hold on
% 
% % Plot concentration profiles
% for idx = 1:length(time_indices)
%     % Red species (warm colors)
%     plot(x, cRed_profiles(time_indices(idx),:),...
%         'LineWidth', 3.5,...
%         'Color', warm_colors(idx,:),...
%         'DisplayName', sprintf('Red t=%.1f s', t(time_indices(idx))));
% 
%     % Ox species (cold colors)
%     plot(x, cOx_profiles(time_indices(idx),:),...
%         'LineWidth', 3.5,...
%         'Color', cold_colors(idx,:),...
%         'DisplayName', sprintf('Ox t=%.1f s', t(time_indices(idx))) );
% end
% 
% xlabel('Distance from Electrode (m)')
% ylabel('Concentration (mmol/L)')
% % title('Concentration Profiles')
% legend('Location', 'northeast')
% % legend off
% grid on
% 
% set(gca,'FontName','Arial','FontSize',50,'FontWeight','bold','LineWidth',2,'Box','off');
% % set(gca,'YScale','log')
%% Core Functions --------------------------------------------------------
function [x,t] = CreateMesh(E, C)
% Creates non-uniform spatial mesh and time vector
% Can be expressed from random walk
L = 60*sqrt(2*abs(E.E_start-E.E_switch)/E.scan_rate*C.D);                     % should this be reactivated??   (YES, ALSO 60 GOOD - In some cases with extreme D's it makes a large difference and improves physicallity of results)  
% L = 0.001;                                                                  % deactivate? -- test! w/ dif Cs
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


% %%
%     function [pl,ql,pr,qr] = bcfun(~,ul,~,ur,t)
%         % Current potential calculation
%         A = E.Electrical.E_switch - E.Electrical.E_start;
%         current_potential = sawtooth(t*pi/abs(A)*...
%             E.Electrical.scan_rate ...
%             , 0.5)*A/2 + (A/2 + E.Electrical.E_start);
% 
%         eta = current_potential - E.Reaction.E0;
% 
% 
% 
%         % Butler-Volmer kinetics
%         alpha = 0.5;
%         k0 = E.Reaction.k0;
%         Fa = E.Reaction.n*C.F*alpha;
%         exp_term = exp(Fa*eta/(C.R*C.T)) ;
%         flux = k0*(ul(1)*exp_term - ul(2)/exp_term);      % Change this to redefine flux
% 
% 
% 
%         % Boundary conditions
%         pl = [flux; flux];  % Proper flux terms for both species
%         ql = [-1; 1];      % Flux directions
%         pr = [0; 0];       % No source terms
%         qr = [1; 1];       % Zero flux conditions
%     end

end

