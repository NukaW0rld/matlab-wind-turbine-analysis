%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ME 4053 PROJECT 2
% PROJECT DESCRIPTION: ANALYSIS OF CLIPPER LIBERTY C96 2.5 MW WIND TURBINE
% AUTHORS: Garrett Kawlewski, James Le, Sophia Steblay
% DATE: 11/01/2024
% DESCRIPTION OF LOCAL VARIABLES:
% 
% %%% Constant Parameters of Wind Turbine
%
% R: Rotor radius [m]
% L: Height of the tower (does not include rotor radius) [m]
% rho: Air density [kg/m^3] at sea level
% A: Rotor swept area [m^2]
% axialFactorMax: Axial induction factor at max power [unitless]
% angularFactorMax: Angular induction factor at max power [unitless]
% mu: Dynamic viscosity of air [kg/(m·s)]
% E: Elastic modulus [Pa]
% sigma_y: Yield strength [Pa]
% S_UT: Ultimate tensile strength [Pa]
% rho_steel: Density of steel [kg/m^3]
%
% %%% Data Imported from Files
%
% blade_data: Table containing blade profile data
% r: Distance from position on blade to axis of rotation [m]
% chord: Chord length at each blade section [m]
% twist: Blade twist angle at each section [deg]
% airfoil_types: Array of airfoil types as strings
% unique_airfoils: Unique airfoil types excluding "circle"
% airfoil_perf: Structure containing airfoil performance data tables
% tower_data: Table containing tower specifications
% h: Height at each tower section [m]
% D_outer: Outer diameter at each tower section [m]
% t: Wall thickness at each tower section [m]
% D_inner: Inner diameter at each tower section [m]
% dP: Power contribution from each blade section [W]
% dT: Thrust contribution from each blade section [N]
% dr: Width of each blade section [m]
% dh: Height of each tower segment [m]
%
% %%% Problem 1 Variables
%
% UProb1: Wind velocity for problem 1 [m/s]
% omegaProb1: Rotational velocity for problem 1 [rad/s]
% pitchProb1: Pitch angle for problem 1 [deg]
% CP_simpson: Coefficient of power for problem 1 [unitless]
% CT_simpson: Coefficient of thrust for problem 1 [unitless]
% PProb1: Calculated power using the coefficient of power [W]
% FTProb1: Calculated thrust force using coefficient of thrust [N]
%
% %%% Problem 2 Variables
%
% UProb2: Wind velocity for problem 2 [m/s]
% tipSpeedRatio2: Tip speed ratio for problem 2 [unitless]
% omegaprob2: Rotational velocity for problem 2 [rad/s]
% beta_sweep_2: Array of pitch angles swept in problem 2 [deg]
% Beta_2: Current pitch angle in loop for problem 2 [deg]
% CP_simpson_2: Coefficient of power for each pitch angle in problem 2 [unitless]
% CT_simpson_2: Coefficient of thrust for each pitch angle in problem 2 [unitless]
% Beta_interp: Interpolated pitch angles for finer resolution [deg]
% CP_interp: Interpolated Coefficient of power [unitless]
% maxValue: Maximum Coefficient of power value [unitless]
% I_max: Index of maximum Coefficient of power [integer]
% Max_pitch_2: Pitch angle corresponding to maximum Coefficient of power [deg]
%
% %%% Problem 3 Variables
%
% UProb3: Wind velocity for problem 3 [m/s]
% pitch_min_Prob3: Minimum pitch angle for problem 3 [deg]
% pitch_max_Prob3: Maximum pitch angle for problem 3 [deg]
% pitch_step_Prob3: Pitch angle step size for problem 3 [deg]
% lambda_min_Prob3: Minimum tip speed ratio for problem 3 [unitless]
% lambda_max_Prob3: Maximum tip speed ratio for problem 3 [unitless]
% lambda_step_Prob3: Tip speed ratio step size for problem 3 [unitless]
% max_CP_Prob3: Maximum Coefficient of Power found in problem 3 [unitless]
% optimal_pitch_Prob3: Optimal pitch angle for maximum coefficient of power [deg]
% optimal_lambda_Prob3: Optimal tip speed ratio for maximum coefficient of power [unitless]
% lambda_values_Prob3: Array of tip speed ratio values for plotting [unitless]
% pitch_angles_Prob3: Array of pitch angle values for plotting [deg]
% CP_matrix_Prob3: Matrix of coefficient of power values over pitch and lambda [unitless]
%
% %%% Problem 4 Variables
% UProb4: Wind velocity for problem 4 [m/s]
% omegaMaxProb4: Maximum rotational velocity for problem 4 [rad/s]
% ratedPower: Rated power of the turbine [W]
% beta_sweep_4: Array of pitch angles swept in problem 4 [deg]
% Beta_4: Current pitch angle in loop for problem 4 [deg]
% CP_simpson_4: Coefficient of Power for each pitch angle in problem 4 [unitless]
% CT_simpson_4: Coefficient of Thrust for each pitch angle in problem 4 [unitless]
% Power_max_pitch: Maximum power achievable without exceeding rated power [W]
% max_blade_pitch_interp: Blade pitch angle corresponding to maximum power [deg]
% Interp_Power: Interpolated power values for finer pitch resolution [W]
% powerDesign: Designed power of the turbine based on Betz limit [W]
% interp_max_power: Maximum interpolated power below rated power [W]
% Pitch_index: Index corresponding to maximum interpolated power [integer]
%
% %%% Problem 5 Variables
% UProb5: Wind speed for problem 5 [m/s]
% omegaProb5: Rotational velocity for problem 5 [rad/s]
% q: Dynamic pressure [N/m^2]
% w: Wind force per unit length at each tower section [N/m]
% ReProb5: Reynolds number at each tower section [unitless]
% CDProb5: Drag coefficient at each tower section [unitless]
% CPProb5: Coefficient of Power for problem 5 [unitless]
% CTProb5: Coefficient of Thrust for problem 5 [unitless]
% FTProb5: Thrust force from the wind turbine [N]
% Torque: Torque generated by the wind turbine [N*m]
% delta_thrust_low: Deflection due to thrust force in low-level model [m]
% delta_wind_low: Deflection due to wind load in low-level model [m]
% delta_static_low: Total static deflection of low-level model [m]
% delta_thrust_high: Cumulative deflection due to thrust in high-level model [m]
% M_wind: Placeholder array for wind moments [N*m]
% beamAMI: Area moment of inertia at each tower section [m^4]
% x_wind_deflection: Array of positions for wind deflection calculation [m]
% beamMod: Modulus of elasticity for the beam [Pa]
% beamL: Length of the beam for deflection calculation [m]
% High_Level_Wind_Deflection: Deflection profile due to wind in high-level model [m]
% Deflection_Wind: Total deflection due to wind at the top of the tower [m]
% delta_static_high: Total static deflection of high-level model [m]
% k_eff_low: Effective stiffness of the low-level model [N/m]
% A_tower: Cross-sectional area of the tower at each section [m^2]
% m_tower_segments: Mass of each tower segment [kg]
% m_tower_total: Total mass of the tower [kg]
% m_eff_low: Effective mass for low-level model using Rayleigh's Method [kg]
% m_nacelle: Mass of the nacelle [kg]
% m_blades: Mass of the blades [kg]
% M_total_low: Total effective mass for low-level model [kg]
% F_test: Test force applied at the top of the tower for stiffness calculation [N]
% v: Deflection array for high-level model [m]
% theta: Slope array for high-level model [rad]
% delta_max_numerical: Maximum deflection from numerical integration in high-level model [m]
% k_eff_high: Effective stiffness of the high-level model [N/m]
% m_eff_high: Effective mass for high-level model using Rayleigh's Method [kg]
% M_total_high: Total effective mass for high-level model [kg]
% omega_n_low: Natural angular frequency for low-level model [rad/s]
% omega_n_high: Natural angular frequency for high-level model [rad/s]
% f_n_low: Natural frequency for low-level model [Hz]
% f_n_high: Natural frequency for high-level model [Hz]
% f0: Fundamental frequency based on rotor speed [Hz]
% pulseWidth: Pulse width for vibration force [s]
% pulseHeight: Pulse height for vibration force [N]
% nHarmonics: Number of harmonics for FFT [unitless]
% fourierTermMag_low: Fourier term magnitudes for low-level model [unitless]
% fourierTermPhase_low: Fourier term phases for low-level model [rad]
% fourierTermMag_high: Fourier term magnitudes for high-level model [unitless]
% fourierTermPhase_high: Fourier term phases for high-level model [rad]
% frequencies: Frequencies corresponding to each harmonic [Hz]
% delta_vibration_low: Dynamic deflection for each harmonic in low-level model [m]
% delta_vibration_high: Dynamic deflection for each harmonic in high-level model [m]
% delta_total_low: Total deflection (static + dynamic) for low-level model [m]
% delta_total_high: Total deflection (static + dynamic) for high-level model [m]
% tower_y: Array of height positions for plotting [m]
% tower_x: Array of x-coordinates for the original tower [m]
% deflection_shape_low_static: Normalized static deflection shape for low-level model [unitless]
% delta_static_low_profile: Static deflection profile for low-level model along height [m]
% delta_vibration_low_profile: Dynamic deflection profile for low-level model along height [m]
% delta_total_low_profile: Total deflection profile for low-level model along height [m]
% deflection_shape_high_thrust: Normalized deflection shape for thrust in high-level model [unitless]
% delta_thrust_high_profile: Deflection profile due to thrust in high-level model along height [m]
% Deflection_Wind_high_interp: Interpolated deflection due to wind for high-level model along height [m]
% delta_vibration_high_profile: Dynamic deflection profile for high-level model along height [m]
% delta_total_high_profile: Total deflection profile for high-level model along height [m]
% Lowest_Inertia: Lowest area moment of inertia for buckling calculation [m^4]
% Min_Area: Minimum cross-sectional area for buckling calculation [m^2]
% Radius_Gyration: Radius of gyration for buckling calculation [m]
% Equivalent_Length_min: Equivalent length for buckling calculation [m]
% Slenderness_Ratio_Part: Slenderness ratio based on equivalent length and radius of gyration [unitless]
% Slenderness_Ratio_JE: Slenderness ratio based on Johnson-Euler buckling [unitless]
% Axial_stress: Axial stress in the tower [MPa]
% Critical_Stress: Critical stress for buckling [MPa]
% Safety_Factor_buckling: Safety factor against buckling [unitless]
% V_x: Shear force at each tower height [N]
% M_thrust_x: Bending moment due to thrust at each tower height [N*m]
% M_wind_x: Bending moment due to wind at each tower height [N*m]
% M_vibration_x: Bending moment due to vibration at each tower height [N*m]
% M_x: Total bending moment at each tower height [N*m]
% Sigma_x: Bending stress at each tower height [MPa]
% Sigma_max: Maximum bending stress [MPa]
% Sigma_max_Static: Maximum bending stress with stress concentration factor [MPa]
% C_g: Surface finish factor [unitless]
% C_r: Reliability factor [unitless]
% C_s: Surface finish factor [unitless]
% S_n: Modified ultimate tensile strength [MPa]
% q: Load factor [unitless]
% K_F: Fatigue factor incorporating stress concentration [unitless]
% Sigma_ea: Effective bending stress for fatigue [MPa]
% Omega_cycles: Rotational speed in cycles [unitless]
% Num_Cycles: Number of cycles the turbine undergoes [cycles]
% S_top: Adjusted ultimate tensile strength for top section [MPa]
% Inf_Life: Infinite life cycle threshold [cycles]
% SF_Fatigue: Safety factor against fatigue failure [unitless]
% Rated_Cycles: Rated number of cycles until failure [cycles]

%
% FUNCTIONS CALLED:
% Calc_Cp_Ct_Main: Calculates the power and thrust coefficients (C_P and C_T) for a wind turbine blade under specified conditions
% cantBeamDistLd: Determine the deflection of a uniformly-loaded cantilever beam
% pulseFFT: Determine the frequency components of a harmonic pulse having constant magnitude
% cylinderCD: To calculate and return the coeffecient of drag for a cyliner in cross flow
% plotSpectrum: Plot frequency content of a square wave pulse
% simpson34: Integrates a function over a given interval using Simpson's 1/3 and 3/8 rules as needed
% sweep_pitch_and_lambda: Finds the optimal blade pitch angle and tip speed ratio for maximum power coefficient at a given wind velocity
% AngleOfAttack_debug: Calculates the relative wind speed (U_rel) and angle of attack (AoA) for a section of the wind turbine blade
%
% START OF EXECUTABLE CODE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTABLISH CONSTANT PARAMETERS OF WIND TURBINE

R = 48; % Rotor radius [m]
L = 77.7; % Height of the tower (does not include rotor radius) [m]
rho = 1.2; % Air density [kg/m^3] at sea level
A = pi * R^2; % Rotor swept area [m^2]
axialFactorMax = 0.33; % Axial induction factor at max power
angularFactorMax = 0.05; % Angular induction factor at max power
mu = 1.8e-5; % Dynamic viscosity of air [kg/(m·s)]
E = 205*10^9; % Elastic modulus [Pa]
sigma_y = 345*10^6; % Yield strength [Pa]
S_UT = 451*10^6; % Ultimate tensile strength [Pa]
rho_steel = 7850; % Density of steel [kg/m^3]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ IN DATA FROM FILES

% Load blade profile data
blade_data = readtable('BladeProfile.csv');
r = blade_data.DistanceFromCenterOfRotation / 1000; % Distance from position on blade to axis of rotation [m]
chord = blade_data.ChordLength / 1000; % Chord length [m]
twist = blade_data.BladeTwist; % Blade twist angle [deg]
airfoil_types = string(blade_data.Airfoil); % Airfoil type as string

% Identify unique airfoil types, excluding the "circle" type that has no corresponding performance data.
unique_airfoils = unique(airfoil_types);
unique_airfoils = unique_airfoils(unique_airfoils ~= "circle");

% Preload airfoil performance data into a structure `airfoil_perf`.
% Note on file naming:
% - The "Airfoil" column in BladeProfile.csv lists names like "DU 97-W-300", "DU 91-W2-250", etc., which contain spaces.
% - However, the actual filenames in the MATLAB directory omit these spaces, e.g., "DU97-W-300" and "DU91-W2-250".
% - This discrepancy is handled by removing spaces from the "Airfoil" entries before constructing the filenames.

airfoil_perf = struct();

for i = 1:length(unique_airfoils)
    % Remove spaces from the airfoil type name to match the actual file names in the directory
    airfoil_name_no_spaces = replace(unique_airfoils(i), ' ', '');
    
    % Convert the modified airfoil name into a valid MATLAB structure field name
    % (using `matlab.lang.makeValidName`) to store the data cleanly within `airfoil_perf`.
    valid_field_name = matlab.lang.makeValidName(unique_airfoils(i));
    
    % Construct the file name by appending ".csv" to the modified airfoil name without spaces
    filename = strcat(airfoil_name_no_spaces, '.csv');
    
    % Check if the file exists in the MATLAB directory before attempting to read it.
    % Note: Originally, these filenames included "V2" at the end (e.g., "DU97-W-300V2.csv"),
    % but the "V2" has been removed from each filename in the directory for compatibility with this script.
    if isfile(char(filename))
        % If the file exists, read its data and store it in the `airfoil_perf` structure
        % under the field name `valid_field_name`.
        airfoil_perf.(valid_field_name) = readtable(char(filename));
    else
        % If the file is not found, throw an error indicating the missing file
        error(['Airfoil performance file ', char(filename), ' not found.']);
    end
end



% Preload tower specs into a structure
tower_data = readtable('towerSpecs.csv');
h = tower_data.Height_mm_ / 1000; % Height [m]
D_outer = tower_data.OD_mm_ / 1000; % Outer diameter [m]
t  = tower_data.WallThk_mm_ / 1000; % Wall thickness [m]
D_inner = D_outer - 2 * t; % Inner diameter [m]

% Initialize arrays for storing section properties
dP = zeros(size(r)); % Power contribution from each section [W]
dT = zeros(size(r)); % Thrust contribution from each section [N]
dr = diff([0; r]); % Width of each section [m]
dh = diff([0; h]); % Segment heights [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM 1: Calculating C_P and C_T for given conditions
% PROBLEM STATEMENT: Determine the coefficient of power and coefficient of
% thrust of the wind turbine for a wind velocity of 10 m/s, rotational
% velocity of 14 rpm, and a pitch angle of 0 degrees.
UProb1 = 10; % Wind velocity [m/s]
omegaProb1 = 14*2*pi/60; % Rotational velocity [rad/s] (converted from rpm)
pitchProb1 = 0; % Pitch angle [rad]

% Call the function with additional output variables
[CP_simpson, CT_simpson] = Calc_Cp_Ct_Main(UProb1, axialFactorMax, angularFactorMax, r, chord, airfoil_types, omegaProb1, pitchProb1, rho, mu,twist);

% Calculate power
PProb1 = CP_simpson * 0.5 * rho * A * UProb1^3;

% Calculate thrust force
FTProb1 = CT_simpson * 0.5 * rho * A * UProb1^2;

% Display the results
fprintf('RESULTS FOR PROBLEM 1:\n');
fprintf('Coefficient of Power (CP_simpson): %.3f\n', CP_simpson);
fprintf('Coefficient of Thrust (CT_simpson): %.3f\n', CT_simpson);
fprintf('Power: %.3f W (%.3f MW)\n', PProb1, PProb1*10^(-6));
fprintf('Thrust force: %.3f N (%.3f kN)\n', FTProb1, FTProb1*10^(-3));
fprintf('Theoretical Maximum C_P (Betz Limit): 0.593\n');
fprintf('Theoretical C_T for a=1/3: 0.889\n');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM 2: Finding max Cp based on pitch angle
% PROBLEM STATEMENT: Sweep the blade pitch angle to determine the maximum
% coefficient of power for wind velocity 9.4 m/s and tip speed ratio 5.88.

fprintf('\nSolving problem 2... (estimated runtime: 5 seconds)\n');

UProb2 = 9.4; % [m/s] Given parameter
tipSpeedRatio2 = 5.88; % [unitless] Given parameter
omegaprob2=(tipSpeedRatio2*UProb2)/(48); %% Finds the rotational velocity by doing the equation (TSR*U/R)

% creates an evenly spaced array of pitch angles from -20 to 30 degrees 
beta_sweep_2=linspace(-20,30,100);
% For loop which creates a pitch array 
for i=1:length(beta_sweep_2)
   Beta_2(i)=beta_sweep_2(i);

[Urel_2, AoA] = AngleOfAttack_debug(UProb2, axialFactorMax, angularFactorMax, r, omegaprob2, Beta_2(i));

% calls the CP and CT calculation equation and stores the return values 
[CP_simpson_2(i), CT_simpson_2(i)] = Calc_Cp_Ct_Main(UProb2, axialFactorMax, angularFactorMax, r, chord, airfoil_types, omegaprob2, Beta_2(i), rho, mu, twist);
%if  CP_simpson_2(i)>0.593
   % CP_simpson_2(i)=0;
% end 

% THIS PRINTS THE VALUE OF CP AND CT AT THE CURRENT ITERATION,ONLY USE FOR DEBUGGING  
% fprintf(' Cp: %f, Ct: %f\n', CP_simpson_2, CT_simpson_2);
end

% Interpolation for more accuracy
% Define new fine pitch angles for interpolation
Beta_interp = linspace(min(Beta_2), max(Beta_2), 100); % 100 points between min and max pitch angles
CP_interp= interp1(Beta_2, CP_simpson_2, Beta_interp, 'linear');

% Find refined maximum value from interpolated data
[maxValue, I_max] = max(CP_interp);
Max_pitch_2 = Beta_interp(I_max); % Corresponding angle for max CP

% Plot 
figure(1);
hold on;
plot(beta_sweep_2,CP_simpson_2);
title('Coeffecient of Power vs Pitch Angle');
xlabel('Pitch Angle (°)');
ylabel('Coeffecient of Power');
xlim([-20 15]);
hold off;
drawnow;
pause(2);

% Output results
fprintf('\nRESULTS FOR PROBLEM 2:\n');
fprintf(' Angle : %.2f, Cp: %.3f\n', Max_pitch_2, maxValue);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM 3: Find max Cp based on pitch angle and tip speed ratio
% PROBLEM STATEMENT: Sweep the blade pitch angle and the tip speed ratio to 
% determine the maximum coefficient of power for wind velocity of 7.6 m/s.

fprintf('\nSolving problem 3... (estimated runtime: 5m30s)\n');

UProb3 = 7.6; % [m/s] Given parameter

% Define pitch angle sweep parameters [degree]
pitch_min_Prob3 = -20;
pitch_max_Prob3 = 30;
pitch_step_Prob3 = 1;

% Define tip speed ratio sweep parameters
lambda_min_Prob3 = 0;
lambda_max_Prob3 = 10;
lambda_step_Prob3 = 0.1;

% Call the sweep_pitch_and_lambda function and stores the max cp found, the
% optimal pitch angle and the optimal tip speed ratio as well as their
% respective arrays for plotting 
[max_CP_Prob3, optimal_pitch_Prob3, optimal_lambda_Prob3, lambda_values_Prob3, pitch_angles_Prob3, CP_matrix_Prob3] = sweep_pitch_and_lambda(UProb3, lambda_min_Prob3, lambda_max_Prob3, lambda_step_Prob3, pitch_min_Prob3, pitch_max_Prob3, pitch_step_Prob3, axialFactorMax, angularFactorMax, r, chord, twist, airfoil_types, airfoil_perf, rho, mu);
pause(1);

% Plots 
figure(2);
hold on;
grid on;
contourf(lambda_values_Prob3, pitch_angles_Prob3, CP_matrix_Prob3, 50, 'LineColor', 'none');
colorbar;
plot(optimal_lambda_Prob3, optimal_pitch_Prob3, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Tip Speed Ratio');
ylabel('Pitch Angle (degrees)');
title(['Maximum CP Sweep for U = ', num2str(UProb3), ' m/s']);
legend('CP Contours', sprintf('Max CP = %.3f at TSR = %.2f, pitch angle = %d°', max_CP_Prob3, optimal_lambda_Prob3, optimal_pitch_Prob3), 'Location', 'best');
hold off;
drawnow;
pause(2);


fprintf('\nRESULTS FOR PROBLEM 3:\n');
fprintf('Max coefficient of power (CP): %.3f\n', max_CP_Prob3);
fprintf('Optimal tip speed ratio: %.2f\n', optimal_lambda_Prob3);
fprintf('Optimal pitch angle: %d degrees\n', optimal_pitch_Prob3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM 4
% PROBLEM STATEMENT: Determine the blade pitch necessary in order to make 
% sure not to over power the turbine between the rated and cut-out speed
% for a wind velocity of 21.2 m/s and a max rotational velocity of 15 rpm

fprintf('\nSolving problem 4... (estimated runtime: 10 seconds)\n');

UProb4 = 21.2; % [m/s] Given parameter
omegaMaxProb4 = 15.5 * 2 * pi() / 60;
ratedPower = 2.5 * 10^6; % Rated power [W]clc 

% evenly spaced array from -10 to 40, range was determined though code
% iteration
 beta_sweep_4=linspace(-20,15,300);

% For loop which creates a pitch array from -10 to 40 degrees 
for i=1:length(beta_sweep_4)
    % New statement with interpolation
     Beta_4(i)=beta_sweep_4(i);

% calls the CP and CT calculation equation and stores the return values 
[CP_simpson_4(i), CT_simpson_4(i)] = Calc_Cp_Ct_Main(UProb4, axialFactorMax, angularFactorMax, r, chord, airfoil_types, omegaprob2,Beta_4(i), rho, mu, twist);

% Conditional statement which sets the Cp to 0 if it exceeds the betz limit
% or is unrealistic
if CP_simpson_4(i)<0 || CP_simpson_4(i)>0.593
   CP_simpson_4(i)=0;
end 
end

% Initialized the value to 0 to avoid syntax errors

Power_max_pitch = 0;


for i=1:1:length(CP_simpson_4)
    
    Power_simpson_4(i) = CP_simpson_4(i) * 0.5 * rho * UProb4^3 * pi() * r(end)^2;

    
    if (Power_simpson_4(i) > Power_max_pitch ) && (Power_simpson_4(i) < ratedPower) && (CP_simpson_4(i) < 1)
            Power_max_pitch = Power_simpson_4(i);
            max_blade_pitch_interp = Beta_4(i);
    end
end
% Finds the power by interpolating
Interp_Power=interp1(Beta_4, Power_simpson_4,beta_sweep_4, 'spline');

% Statement which locates the highest power thats less than the rated power
interp_max_power = max(Interp_Power(Interp_Power < ratedPower));
Pitch_index=find(Interp_Power == interp_max_power, 1);
max_blade_pitch_interp = beta_sweep_4(Pitch_index);

% Calculate CT using U and pitch angle
[CPProb5, CTProb5] = Calc_Cp_Ct_Main(UProb5, axialFactorMax, angularFactorMax, r, chord, airfoil_types, omegaProb5, max_blade_pitch_interp, rho, mu, twist);

FTProb5 = CTProb5 * 0.5 * rho * UProb5^2 * A; % Thrust force [N]



fprintf('\nRESULTS FOR PROBLEM 4:\n');
fprintf('Blade Pitch Angle [degrees]: %.2f, Power at Blade Pitch Angle [MW]: %.2f\n', max_blade_pitch_interp, interp_max_power * 10^-6);
fprintf('Thrust force: %.3f N\n', FTProb5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROBLEM 5

UProb5 = 21.2; % Wind speed [m/s]
omegaProb5 = 15.5 * 2*pi/60; % Rotational velocity [rad/s]

% 1. STATIC WIND LOADS ON THE TOWER
q = 0.5 * rho * UProb5^2; % Dynamic pressure [N/m^2]

w = zeros(length(h), 1); % Wind force per unit length [N/m]
ReProb5 = zeros(length(h), 1); % Reynolds number
CDProb5 = zeros(length(h), 1); % Drag coefficient

% Calculate w(x) for each height segment [wind force per unit length; N/m]
for i = 1:length(h)
    ReProb5(i) = (rho * UProb5 * D_outer(i)) / mu; % Reynolds number
    CDProb5(i) = cylinderCD(ReProb5(i)); % Calculate CD from Reynolds number
    w(i) = q * CDProb5(i) * D_outer(i); % Wind force per unit length
end

% 2. DETERMINE TORQUE


Torque=CPProb5* 0.5 * rho * UProb5^3 * pi * r(end)^2/omegaProb5;


% 3. DO STATIC DEFLECTION ANALYSIS

% Low-level model
D_outer_avg = mean(D_outer);
t_avg = mean(t);
I_avg = pi/64 * (D_outer_avg^4 - (D_outer_avg - 2*t_avg)^4); % Average area moment of inertia [m^4]
w_avg = mean(w); % Average distributed load [N]

delta_thrust_low = (FTProb5 * L^3)/(3 * E * I_avg); % Deflection due to thrust force [m]
delta_wind_low = (w_avg * L^4)/(8 * E * I_avg); % Deflection due to wind [m]


delta_static_low = delta_thrust_low + delta_wind_low; % Total deflection of the low-level model [m]

% High-level model
delta_thrust_high = 0;
% delta_wind_high = 0;
M_wind = zeros(length(h), 1);
beamAMI = zeros(length(h), 1);
x_wind_lowlevel=linspace(0,h(end),length(h));



for i=1:length(h)
    beamAMI(i) = pi/64 * (D_outer(i)^4 - D_inner(i)^4); % Area moment of inertia at x [m^4]
end 

x_wind_deflection=linspace(0,h(end),length(beamAMI)); % Linearly spaced array
beamMod=E; %Pa 


% w=q; % N/m
beamL=h(end); % m

% assigns the return value of cantBeamDistLd to y_p1
High_Level_Wind_Deflection=cantBeamDistLd(x_wind_deflection,w,beamL,beamMod,beamAMI);
Deflection_Wind=High_Level_Wind_Deflection(end);

for i = 1:length(h) % Calculate deflection using numerical integration
    x = h(i); % Current position [m]
    I_x = pi/64 * (D_outer(i)^4 - D_inner(i)^4); % Area moment of inertia at x [m^4]

    % Deflection due to thrust force [m]
    delta_thrust_high = delta_thrust_high + (FTProb5 * (L - x)^2) / (E * I_x) * dh(i);
end

delta_static_high = delta_thrust_high + Deflection_Wind; % Total deflection of the high-level model [m]

fprintf('\nRESULTS FOR PROBLEM 5\n');
fprintf('The static deflection of the low-level model is: %.6f m\n', delta_static_low);
fprintf('The static deflection of the high-level model is: %.6f m\n', delta_static_high(end));

%% 4. DYNAMIC VIBRATION ANALYSIS

% Effective mass and stiffness calculation
% Low-level model

% Calculate effective stiffness
k_eff_low = (3 * E * I_avg) / L^3; % [N/m]

% Calculate total mass of the tower
% Mass per segment: rho_steel * Area * height
A_tower = pi * (D_outer.^2 - D_inner.^2) / 4; % cross-sectional area [m^2]
m_tower_segments = rho_steel .* A_tower .* dh; % mass per segment [kg]
m_tower_total = sum(m_tower_segments); % total mass [kg]

% Calculate effective mass (m_eff_low) using Rayleigh's Method for low-Level model
% For a cantilever beam with uniform mode shape, m_eff = m_tower_total / 3
m_eff_low = m_tower_total / 3; % [kg]

% Calculate total effective mass (m_total_low) for low-level model
m_nacelle = 94.3e3; % mass of nacelle [kg]
m_blades = 55.7e3; % mass of blades [kg] (all 3 blades already included; do not multiply by 3)
M_total_low = m_nacelle + m_blades + m_eff_low; % total effective mass [kg]

% High-level model

% Calculate effective stiffness (k_eff_high) for high-Level model using numerical method

% Define a test force (F_test) to apply at the top of the tower
F_test = 1; % [N] (arbitrary small force for stiffness calculation)

% Initialize deflection and slope arrays
v = zeros(length(h)+1, 1); % deflection [m] (length+1 for boundary condition)
theta = zeros(length(h)+1, 1); % slope [rad]

% Numerical integration using Euler's Method
for i = 1:length(h)
    x = h(i); % current height [m]
    % Bending moment at position x due to F_test
    M_i = F_test * (L - x); % [N*m]
    % Curvature (d^2v/dx^2) = -M / (E * I(x))
    curvature = -M_i / (E * (pi/64) * (D_outer(i)^4 - D_inner(i)^4)); % [1/m]
    % Update slope (theta) using curvature
    theta(i+1) = theta(i) + curvature * dh(i); % [rad]
    % Update deflection (v) using slope
    v(i+1) = v(i) + theta(i) * dh(i); % [m]
end

% Maximum deflection at the top (delta_max_numerical)
delta_max_numerical = v(end); % [m]

% Apply absolute value to ensure positive deflection
delta_max_numerical = abs(delta_max_numerical);

% Calculate effective stiffness (k_eff_high)
k_eff_high = F_test / delta_max_numerical; % [N/m]

% Calculate effective mass (m_eff_high) using Rayleigh's Method for high-Level model

% Define mode shape for first bending mode of a cantilever beam
% v(x) = (x/L)^2 * (3 - 2*x/L) * v_max
% Normalize mode shape by v_max to make it dimensionless
v_ratio = (h/L).^2 .* (3 - 2*(h/L)); % dimensionless mode shape

% Calculate effective mass using Rayleigh's Method
m_eff_high = sum(rho_steel * A_tower .* (v_ratio.^2) .* dh); % [kg]

% Calculate total effective mass (M_total_high) for high-Level model
M_total_high = m_nacelle + m_blades + m_eff_high; % total effective mass [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Natural frequency calculation

% Low-level model

% Calculate natural angular frequency (omega_n) for low-Level model
omega_n_low = sqrt(k_eff_low / M_total_low); % [rad/s]
omega_n_high = sqrt(k_eff_high / M_total_high); % [rad/s]

% Convert to natural frequency (f_n)
f_n_low = omega_n_low / (2 * pi); % [Hz]
f_n_high = omega_n_high / (2 * pi); % [Hz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harmonic decomposition of vibration force

% Define vibration force characteristics
% Assuming vibration force is periodic with fundamental frequency equal to rotor frequency
f0 = omegaProb5 / (2 * pi); % fundamental frequency [Hz]

% Define pulse characteristics for vibration force
pulseWidth = 0.09 / f0; % pulse width as a fraction of period (must be less than 0.1/f0)
pulseHeight = FTProb5; % approximate vibration force magnitude [N]
nHarmonics = 100; % number of harmonics to include (must be less than 500)

% Generate harmonic components using pulseFFT.m
[fourierTermMag_low, fourierTermPhase_low] = pulseFFT(f0, pulseWidth, pulseHeight, nHarmonics);
[fourierTermMag_high, fourierTermPhase_high] = pulseFFT(f0, pulseWidth, pulseHeight, nHarmonics);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency spectrum visualization

% Plot deflection response spectrum for low-level model
plotSpectrum(f0, pulseWidth, nHarmonics, fourierTermMag_low, 3);

% Plot deflection response spectrum for high-level model
plotSpectrum(f0, pulseWidth, nHarmonics, fourierTermMag_high, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displacement response calculation

% Frequencies corresponding to each harmonic
frequencies = (1:nHarmonics) * f0; % [Hz]

% Initialize displacement arrays
delta_vibration_low = zeros(nHarmonics, 1); % dynamic deflection for each harmonic [m]
delta_vibration_high = zeros(nHarmonics, 1); % dynamic deflection for each harmonic [m]

for n = 1:nHarmonics
    % Angular frequency for nth harmonic
    omega_n_harmonic = 2 * pi * n * f0; % [rad/s]
    
    % Low-level model response (undamped)
    denominator_low = abs(k_eff_low - M_total_low * omega_n_harmonic^2); % [N/m]
    delta_vibration_low(n) = fourierTermMag_low(n) / denominator_low;
    
    % High-level model response (undamped)
    denominator_high = abs(k_eff_high - M_total_high * omega_n_harmonic^2); % [N/m]
    delta_vibration_high(n) = fourierTermMag_high(n) / denominator_high;
end

% Plot deflection response spectrum for low-level model
figure(5);
stem(frequencies, delta_vibration_low, 'b');
xlabel('Frequency (Hz)');
ylabel('Deflection amplitude (m)');
title('Deflection response spectrum - Low-level model');
grid on;

% Plot deflection response spectrum for high-level model
figure(6);
stem(frequencies, delta_vibration_high, 'r');
xlabel('Frequency (Hz)');
ylabel('Deflection amplitude (m)');
title('Deflection response spectrum - High-level model');
grid on;

% Display natural frequencies
fprintf('\nLow-level model natural frequency: %.2f Hz\n', f_n_low);
fprintf('High-level model natural frequency: %.2f Hz\n', f_n_high);

% Display dynamic deflection results
% Report maximum deflection response
delta_vibration_max_low = max(abs(delta_vibration_low));
delta_vibration_max_high = max(abs(delta_vibration_high));

fprintf('Maximum deflection due to vibration (low-level model): %.6f m\n', delta_vibration_max_low);
fprintf('Maximum deflection due to vibration (high-level model): %.6f m\n', delta_vibration_max_high);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resonance assessment

% Define a threshold for resonance (within 5% of natural frequency)
threshold = 0.05; % 5%

% Low-level model resonance check
resonance_low = false;
for n = 1:nHarmonics
    harmonic_freq = n * f0;
    if f_n_low ~= 0 && abs(harmonic_freq - f_n_low) / f_n_low < threshold
        fprintf('Resonance alert: harmonic %d at %.2f Hz is within %.2f%% of natural frequency (low-level model).\n', ...
            n, harmonic_freq, threshold*100);
        resonance_low = true;
    end
end
if ~resonance_low
    disp('No resonance detected in low-level model.');
end

% High-level model resonance check
resonance_high = false;
for n = 1:nHarmonics
    harmonic_freq = n * f0;
    if f_n_high ~= 0 && abs(harmonic_freq - f_n_high) / f_n_high < threshold
        fprintf('Resonance alert: harmonic %d at %.2f Hz is within %.2f%% of natural frequency (high-level model).\n', ...
            n, harmonic_freq, threshold*100);
        resonance_high = true;
    end
end
if ~resonance_high
    disp('No resonance detected in high-level model.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total deflection calculation

% Total deflection including vibration
delta_total_low = delta_static_low + delta_vibration_max_low;
delta_total_high = delta_static_high + delta_vibration_max_high;

fprintf('\nTotal deflection (low-level model): %.6f m\n', delta_total_low);
fprintf('Total deflection (high-level model): %.6f m\n', delta_total_high);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting total deflection


tower_y = linspace(0, L, 50); % Make a plot of the tower without being deflected, going straight up from x = 0 to y = L
tower_x = zeros(size(tower_y));

% Low-level model deflection

% For a cantilever beam under point load at the tip, deflection shape is:
% delta(y) = delta_static_low * (y^2(3L - y))/L^3
deflection_shape_low_static = (tower_y.^2 .* (3*L - tower_y)) / (L^3);
delta_static_low_profile = delta_static_low * deflection_shape_low_static; % [m]
delta_vibration_low_profile = delta_vibration_max_low * (tower_y / L); % [m]
delta_total_low_profile = delta_static_low_profile + delta_vibration_low_profile;

% High-level model deflection

% Since delta_thrust_high is a scalar representing the total deflection due
% to thrust at the top, we have to calculate deflection due to thrust at
% each y by distributing the thrust deflection similar to low-level model
deflection_shape_high_thrust = (tower_y.^2 .* (3*L - tower_y)) / (L^3);
delta_thrust_high_profile = delta_thrust_high * deflection_shape_high_thrust; % [m]

% We'll interpolate Deflection_Wind to match tower_y
Deflection_Wind_high_interp = interp1(x_wind_deflection, High_Level_Wind_Deflection, tower_y, 'linear', 'extrap');

delta_vibration_high_profile = delta_vibration_max_high * (tower_y / L); % [m]

% Total deflection for high-level model
delta_total_high_profile = Deflection_Wind_high_interp + delta_thrust_high_profile + delta_vibration_high_profile;

% Plotting
figure(7);
hold on;
grid on;
plot(tower_x, tower_y, 'k-', 'LineWidth', 2, 'DisplayName', 'Original tower'); % Plot the original undeflected tower
plot(delta_total_low_profile, tower_y, 'LineWidth', 2, 'DisplayName', 'Low-level model'); % Plot the low-level model
plot(delta_total_high_profile, tower_y, 'LineWidth', 2, 'DisplayName', 'High-level model'); % Plot the high-level model
xlabel('Deflection (m)');
ylabel('Height (m)');
title('Wind Turbine Tower Deflection Under Loading Condition U = 21.2 m/s');
legend('Location', 'best');

xlim([-3, 3]);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating vibration force
F_eff_low = k_eff_low * abs(delta_vibration_max_low); % [N]
F_eff_high = k_eff_high * abs(delta_vibration_max_high); % [N]







%% 5. DO STRESS ANALYSIS TO EVALUATE STRENGTH
% Buckling
% Calculates the inertia for each iteration and finds the
% lowest inertia for buckling calculation
    Lowest_Inertia=(D_outer(end)^4-D_inner(end)^4)*pi/64;
    Min_Area=(pi*(D_outer(end)^2-D_inner(end)^2))/4;

  
Radius_Gyration=sqrt(Lowest_Inertia/Min_Area);
% Assumed Fixed-Open Beam, like Cantilever beam 
Equivalent_Length_min=2*h(end);
 % Le/rho of the part calculation 
 Slenderness_Ratio_Part=Equivalent_Length_min/Radius_Gyration;
 % Le/rho of Johnson-Euler calculation 
 Slenderness_Ratio_JE=sqrt((2*pi^2*E)/sigma_y);
 Axial_stress=(1.4*(m_nacelle+m_blades)*9.81)/(Min_Area*1E6);

 % Condition which determines if Johnson or Euler is needed 
 if Slenderness_Ratio_JE>Slenderness_Ratio_Part
     Critical_Stress=abs(((sigma_y-(sigma_y^2*Slenderness_Ratio_Part^2))/(4*pi^2*E)));
     % Need the force here acting on the y plane i think
     
     Safety_Factor_buckling=(Critical_Stress*Min_Area)/((m_nacelle+m_blades)*9.81);
 else 
     Critical_Stress=(pi^2*E)/Slenderness_Ratio_Part^2;
    Safety_Factor_buckling=Critical_Stress*Min_Area/((m_nacelle+m_blades)*9.81);
 end 
   fprintf('\n The Axial Stress is: %.2f MPa\n', Axial_stress);

 fprintf('\n The Safety Factor of buckling on the smallest cross-section is: %.2f\n', Safety_Factor_buckling);



%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stress Analysis
% Initialization of each array 
V_x = zeros(length(h), 1);
M_thrust_x = zeros(length(h), 1);
M_wind_x = zeros(length(h), 1);
M_vibration_x = zeros(length(h), 1);
M_x = zeros(length(h), 1);
sigma_bend_x = zeros(length(h), 1);
Kt = 3; % Static Stress concentration factor due to flanges,


for i = 1:length(h)
    x = h(i); % Current height [m]
    s = h(i:end); % Positions from current x to L [m]
    c = D_outer(i) / 2; % Distance from neutral axis [m]

    % Shear force caused by thrust force and wind force
    V_x(i) = FTProb5 + sum(w(i:end) .* dh(i:end)); % [N]

    % Bending moments caused by thrust force, wind force, and vibration force
    M_thrust_x(i) = FTProb5 * (L - x); % Moment due to thrust force acting at the top of the tower [N*M]
    M_wind_x(i) = sum(w(i:end) .* dh(i:end) .* (s - x)); % Moment due to wind force acting as distributed load along the tower [N*m]
    M_vibration_x(i) = F_eff_high * (L - x); % [N*m]
    
    % Total bending moment at height x
    M_x(i) = M_thrust_x(i) + M_wind_x(i) + M_vibration_x(i);

    % Area moment of inertia
    Inertia=pi*(D_outer(i)^4-D_inner(i)^4)/64;

    % Bending stress caused by bending moments
    Sigma_x(i)=(M_x(i)*c/(Inertia*1E6)); % [MPa]
    % fprintf('\n The bending stress is: %.2f Mpa\n', Sigma_x_total);
end

% Plot the bending stress
figure(8);
hold on;
grid on;
plot(h,Sigma_x*Kt)
xlabel('Tower height (m)');
ylabel('Amplified bending stress (MPa)');
title('Bending stress amplified by stress concentration factor Kt = 3 along the height of the tower');
hold off;

% Calculate the bending stress, the static stress and outputs the axial
% stress, the safety factor of the material, and the ratio between axial and 
Sigma_max=max(Sigma_x);
Sigma_max_Static=Sigma_max*Kt;
fprintf('\n The bending stress is: %.2f Mpa\n', Sigma_max_Static);
fprintf('\n The Safety Factor for static bending stress is: %.2f\n', (sigma_y/(Sigma_max_Static*1E6)));
fprintf('\n The Axial stress is: %.2f Percent of the bending stress \n', 100*Axial_stress/Sigma_max_Static);


%% Fatigue Analysis 
C_g=0.6; % Assuming diameter>250mm
C_r=0.814; % Reliability factor, assuming 99% reliability 
C_s = 0.8; % Assuming surface finish factor for machined surfaces
S_n=C_g*C_r*C_s*0.5*S_UT*1E-6; % Endurance strength, assuming machined finished, 99% reliability, and and a diameter >250 mm, in MPa
q=.75; % Notch Sensitivity
K_F=1+q*(Kt-1); % Dynamic stress concentration factor calculation, needed for fully reversible loading 
Sigma_ea=(Sigma_max)*K_F; % Alternating stress calculation, assuming 
Omega_cycles=15.5;
% Calculates cycles 
Num_Cycles=15.5*60*24*365*3; % Number of cycles for 3 years [Cycles/3 years]
% Fully reversed loading
S_top=0.9*S_UT; % 10^3 Cycle strength
% Infinite life cycles for steel
Inf_Life=10^6;
% Basquin's constant for steel
b=-.1;




%Conditional state which calculates the safety factor for fully reversible
%loading 
if Num_Cycles>Inf_Life
SF_Fatigue=(S_n/Sigma_ea);

% If the 
if SF_Fatigue<1
  fprintf('\n The Safety Factor is: %.2f ,Infinite life is not achievable \n', SF_Fatigue);

  Rated_Cycles=(SF_Fatigue)^(1/b);
   fprintf('\n The number of cycles until failure is : %2d \n', Rated_Cycles);

end
end

















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FUNCTIONS USED IN THE PROGRAM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









function [CP_simpson, CT_simpson] = Calc_Cp_Ct_Main(U, alpha, alpha_p, r, chord, airfoil_types, omega, pitch, rho, mu, twist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:
%   Calc_Cp_Ct_Main
%
% PURPOSE:
%   Calculates the Coefficient of Power (C_P) and Coefficient of Thrust (C_T)
%   for a wind turbine blade under specified conditions. The function 
%   evaluates power and thrust contributions along the blade span based on 
%   local blade properties and environmental parameters.
%
% INPUT:
%   U - Wind velocity [m/s]
%   alpha - Axial induction factor [unitless]
%   alpha_p - Angular induction factor [unitless]
%   r - Radial positions along the blade [m]
%   chord - Chord lengths at local radial position [m]
%   airfoil_types - Array of airfoil types corresponding to radial positions
%   omega - Rotational velocity of the blade [rad/s]
%   pitch - Blade pitch angle [deg]
%   rho - Air density [kg/m^3]
%   mu - Dynamic viscosity of air [kg/(m·s)]
%   twist - Blade twist angle [deg]
%
% OUTPUT:
%   CP_simpson - Coefficient of Power
%   CT_simpson - Coefficient of Thrust
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure pitch is a vector of the correct size
if isscalar(pitch)
    pitch = repmat(pitch, 1, length(r));
end 

% Initialize arrays
CL = zeros(1, length(r));
CD = zeros(1, length(r));
Fz_array = zeros(1, length(r));
F_theta_array = zeros(1, length(r));
DU97 = readtable('DU97-W-300.csv', 'NumHeaderLines', 1);
DU91 = readtable('DU91-W2-250.csv', 'NumHeaderLines', 1);
DU93 = readtable('DU93-W-210.csv', 'NumHeaderLines', 1);
DU96 = readtable('DU96-W-180.csv', 'NumHeaderLines', 1);

% Loop through each blade section
for i = 1:length(r)
    local_r = r(i);
    local_chord = chord(i);
    theta_P = deg2rad(pitch(i)) + deg2rad(twist(i)); 
    U_axial = U * (1 - alpha);          
    U_radial = local_r * omega * (1 - alpha_p);  
    U_rel = sqrt(U_axial^2 + U_radial^2);
    phi = atan(U_axial / U_radial);
    AoA = phi - theta_P;


    % Load airfoil performance data based on radial position
   
        if airfoil_types(i)== "circle"
            CL(i) = 0;
            Re = rho * U_rel * local_chord / mu;
            CD(i) = cylinderCD(Re); % Example drag value for non-lift region
        elseif airfoil_types(i)== "DU 97-W-300"
           % DU97 = readtable('DU97-W-300.csv', 'NumHeaderLines', 1);
            AoA_data = deg2rad(table2array(DU97(:, 1))');
            cL_data = table2array(DU97(:, 2))';
            cD_data = table2array(DU97(:, 3))';
            CL(i) = interp1(AoA_data, cL_data, AoA, 'linear', 'extrap');
            CD(i) = interp1(AoA_data, cD_data, AoA, 'linear', 'extrap');
        elseif airfoil_types(i)== "DU 91-W2-250"
            %DU91 = readtable('DU91-W2-250.csv', 'NumHeaderLines', 1);
            AoA_data = deg2rad(table2array(DU91(:, 1))');
            cL_data = table2array(DU91(:, 2))';
            cD_data = table2array(DU91(:, 3))';
            CL(i) = interp1(AoA_data, cL_data, AoA, 'linear', 'extrap');
            CD(i) = interp1(AoA_data, cD_data, AoA, 'linear', 'extrap');
        elseif airfoil_types(i)=="DU 93-W-210"
            % DU93 = readtable('DU93-W-210.csv', 'NumHeaderLines', 1);
            AoA_data = deg2rad(table2array(DU93(:, 1))');
            cL_data = table2array(DU93(:, 2))';
            cD_data = table2array(DU93(:, 3))';
            CL(i) = interp1(AoA_data, cL_data, AoA, 'linear', 'extrap');
            CD(i) = interp1(AoA_data, cD_data, AoA, 'linear', 'extrap');
       elseif airfoil_types(i)== "DU 96-W-180"
          %  DU96 = readtable('DU96-W-180.csv', 'NumHeaderLines', 1);
            AoA_data = deg2rad(table2array(DU96(:, 1))');
            cL_data = table2array(DU96(:, 2))';
            cD_data = table2array(DU96(:, 3))';
            CL(i) = interp1(AoA_data, cL_data, AoA, 'linear', 'extrap');
            CD(i) = interp1(AoA_data, cD_data, AoA, 'linear', 'extrap');
        end

    % Calculate forces
    dL = CL(i) * 0.5 * rho * U_rel^2 * local_chord;
    dD = CD(i) * 0.5 * rho * U_rel^2 * local_chord;
    psi = atan((1 - alpha) * U / ((1 - alpha_p) * omega * local_r));

    Fz_array(i) = dL * cos(psi) + dD * sin(psi); % Axial force
    F_theta_array(i) = (dL * sin(psi) - dD * cos(psi)) * local_r; % Tangential force
end

% Calculate available wind power and thrust
wind_power = 0.5 * rho * U^3 * pi * r(end)^2;
wind_thrust = 0.5 * rho * U^2 * pi * r(end)^2;

% Interpolation for Simpson's rule
r_new = linspace(r(1), r(end), 100);
Torquet_interp = interp1(r, F_theta_array, r_new, 'spline');
Forcen_interp = interp1(r, Fz_array, r_new, 'spline');

% Calculate torque and thrust using Simpson's rule
Torquet_simpson = 3 * simpson34(r_new(1), r_new(end), length(r_new), r_new, Torquet_interp);
Thrustn_simpson = 3 * simpson34(min(r_new), max(r_new), length(r_new), r_new, Forcen_interp);

% Total power and thrust calculations
total_power_simpson = Torquet_simpson * omega;
total_thrust_simpson = Thrustn_simpson;

% Coefficients of Power and Thrust
CP_simpson = total_power_simpson / wind_power;
CT_simpson = total_thrust_simpson / wind_thrust;

end





function y = cantBeamDistLd(x, w, beamL, beamMod, beamAMI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME: cantBeamDistLd
%
% PURPOSE
% Determine the deflection of a uniformly-loaded cantilever beam
%
% INPUT
% x - Array of x-positions along beam where deflections should be returned
% w - Intensity of applied load (N/m)
% beamL - Length of beam (meters)
% beamMod - Modulus of elasticity (Pa)
% beamAMI - Area moment of inertia of beam cross-section (m^4)
%             
%
% OUTPUT
% y - Array of beam deflections at input positions x (meters)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Ensure x is a column vector
    x = x(:);
    
    if isscalar(beamAMI)     
        I = @(x_val) beamAMI; % If beamAMI is a scalar, define I(x) as a constant function
    else
        I = @(x_val) interp1(linspace(0, beamL, length(beamAMI)), beamAMI, x_val, 'linear', 'extrap');
        % If beamAMI is an array, create an interpolation function
    end

    % ((C)) Handle variable load w
    if isscalar(w)
        w_func = @(x_val) w; % Constant load
    else
        w_func = @(x_val) interp1(linspace(0, beamL, length(w)), w, x_val, 'linear', 'extrap');
    end
    % End of ((C))
    
    % Initial conditions at x = 0 (fixed end)
    y0 = [0; 0]; 

    odefun = @(x_val, y_vec) [y_vec(2); -w_func(x_val) * (beamL * x_val - (x_val.^2) / 2) / (beamMod * I(x_val))]; 
    [xsol, ysol] = ode45(odefun, [0, max(x)], y0);
    y = interp1(xsol, ysol(:,1), x, 'linear', 'extrap');    
end






function [ fourierTermMag, fourierTermPhase  ] = ...
   pulseFFT( freq, pulseWidth, pulseHeight, nHarmonics )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: pulseFFT
%
%  PURPOSE
%    Determine the frequency components of a harmonic pulse
%    having constant magnitude
%
%  INPUT
%    freq: Frequency of the input pulses (Hz)
%    nHarmonics: # of harmonics to include in the Fourier series
%                NOTE 1: The fundamental frequency is counted as
%                        the "first" harmonic
%                NOTE 2: nHarmonics must be < 500 due to the
%                        sampling rate assumed in this function
%    pulseHeight: Magnitude of the pulse 
%                 NOTE: Units must be consistent with forceMag
%    pulseWidth: Width of the pulse (sec)
%                NOTE: pulseWidth must be less than 1/freq!
%
%  OUTPUT
%    fourierTermMag: Magnitudes of the harmonic terms 
%                    (same units as pulseHeight)
%            NOTE 1: Units must be consistent with pulseHeight
%            NOTE 2: forceMag(1) = DC component
%                    forceMag(2) = magnitude of fundamental freq
%                    forceMag(3) = magnitude of 1st harmonic
%                                  (2X fundamental freq)
%                    forceMag(i) = magnitude of (i-2)nd harmonic
%                                  ((i-2)X fundamental freq)
%   fourierTermPhase: Phase angles of the harmonic terms (radians)
%                     (see Note 2 above for meaning of each term)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: F. Kelso & T. Chase
%  DATE: 10/25/24
%
%  DESCRIPTION OF LOCAL VARIABLES
%    ampPhase: Complex number representation of the raw
%              amplitude and phase of one harmonic component
%              of the input pulse
%    denom: Scaling factor to convert the raw digital representation
%           of the discrete Fourier transform of the input signal
%           to a usable magnitude and phase angle
%    index: Loop counter index
%    nPts: Number of points used to digitally sample the pulse
%    nPulse: Number of points of the digital signal 
%            that the pulse is "on" (between 0 and nPts)
%    pulseFrac: Fraction of the signal period that the pulse is "on"
%    samplePeriod: Period of the digital sample
%    signal: Raw discrete Fourier transform of the input pulse
%    sigPeriod: Period of the input signal (s)
%    sqrPulse: Array representing the pulse as a function of time
%              over the period of the input signal
%
%  FUNCTIONS CALLED
%    fft: Standard MATLAB function for performing a
%         discrete Fourier transform of an input vector
%    floor: Standard MATLAB function for rounding down
%           to the nearest integer toward infinity of a real number
%    round: standard MATLAB function for rounding a real number
%           to the closest integer
%    zeros: standard MATLAB function for creating an array of zeros
%
%  START OF EXECUTABLE CODE
%
%  Check for a legal number of harmonics
%
   if nHarmonics >= 500
%
%  Too many harmonics requested.
%  All return arguments must be set to something to avoid a MATLAB error.
%
      fourierTermMag = NaN;
      fourierTermPhase = NaN;
      return
   end
%
%  Determine period of cyclical pulse
%
   sigPeriod = 1.0 / freq;
%
%  Check for a legal pulse width
%
   if pulseWidth >= sigPeriod
%
%  Input is not a pulse!
%  All return arguments must be set to something to avoid a MATLAB error.
%
      fourierTermMag = NaN;
      fourierTermPhase = NaN;
      return
   end
%
%  Assume that signal is sampled using 1000 data points
%
   nPts = 1000;  
%
%  Determine period of sample
%
   samplePeriod = sigPeriod / nPts;
%
%  Determine fraction of the period where pulse is "on"
%
   pulseFrac = pulseWidth / sigPeriod; 
%
%  Define an array defining the pulse as a function of time...
%
%  ...Start with "off"
%
   sqrPulse = zeros( 1, nPts );
%
%  ...Now add the "on"
%     (+1 accounts for the pulse being "on" at zero time)
%
   nPulse = round( pulseFrac * nPts ) + 1; 
   sqrPulse( 1:nPulse ) = pulseHeight;
%
%  Perform a discrete Fourier transform on the input pulse

   signal = fft( sqrPulse ); %(complex) components of the input pulse
%
%  Convert the raw complex components of the input pulse
%  to magnitudes and angles...
%  (A reference needs to be added here to explain the theory behind
%   what is being done)
%
%  ...Reserve space for the output arrays
%
   fourierTermMag = zeros( 1, nHarmonics );
   fourierTermPhase = zeros( 1, nHarmonics );
%
%  ...Determine the scale factor
%
   denom = floor( nPts );
%
%  ...Set the magnitude and phase angle (0) of the DC component
%     of the signal
%
   fourierTermMag( 1 ) = signal( 1, 1 ) / denom;
   fourierTermPhase( 1 ) = 0.0;
%
%  ...Calculate all harmonics requested by the user
%
   for index = 2:nHarmonics
%
%  ...Capture the complex # representation of one harmonic
%
      ampPhase = 2.0 * signal( index ) / denom;
%
%  ...Convert to magnitude and phase angle
%
      fourierTermMag( index ) = abs( ampPhase );
      fourierTermPhase( index ) = angle( ampPhase );
   end
end







function [C_D] = cylinderCD(Re)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: cylinderCD
%
%  PURPOSE:
%   To calculate and return the coeffecient of drag for a cyliner in cross
%   flow.  
%
%   Fit functions can be found in references:
%   
%   - White, F.M. (2006) Viscous Fluid Flow. 3rd Edition, McGraw-Hill, Boston.
%   - Nian-Sheng Cheng, *Calculation of Drag Coefficient for Arrays of ...
%      Emergent Circular Cylinders with Pseudofluid Model.*...
%      J. Hydraul. Eng. 2013.139:602-611.
%  
% INPUT
%   Re - Reynolds number [-]
%
%  OUTPUT
%   C_D - drag coeffecient [-]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: MJA
%  DATE: 2024.10.21
%
%  DESCRIPTION OF LOCAL VARIABLES
%  
%
%  FUNCTIONS CALLED
%   log10 (MATLAB)
%   tanh (MATLAB)
%
%  START OF EXECUTABLE CODE
%

% piecewise fit functions based on various reynolds number regimes for a
% cylinder in cross-flow

if Re < 2*10^5
    C_D = 11 * Re.^(-0.75) + 0.9 * (1.0 - exp(-1000./Re))...
        + 1.2 * (1.0 -exp(-(Re./4500).^0.7));

elseif Re <= 5*10^5
    C_D = 10.^(0.32*tanh(44.4504 - 8 * log10 (Re))-0.238793158);

else
    C_D = 0.1 * log10(Re) - 0.2533429;

end
end









function plotSpectrum( freq, pulseWidth, nHarmonics, harmMag, figNo )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: plotSpectrum
%
%  PURPOSE
%    Plot frequency content of a square wave pulse
%
%  INPUT
%    freq: Frequency of the cyclic signal (Hz)
%    pulseWidth: Width of the square pulse (s)
%    nHarmonics: Number of harmonics used to represent the pulse
%                (includes DC & fundamental frequency)
%    harmMag: Array of magnitudes of pulses at each harmonic
%    figNo: Number to be assigned to plot
%
%  OUTPUT
%    None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: TRC
%  DATE: 10/26/24
%
%  DESCRIPTION OF LOCAL VARIABLES
%
%  FUNCTIONS CALLED
%    axis: Standard MATLAB function to contol the scaling and appearance
%          of a plot
%    figure: Standard MATLAB function to create a figure window
%    num2str: Standard MATLAB function to convert numbers to characters
%    stem: Standard MATLAB function to create a discrete sequence
%          or "stem" plot
%    subtitle: Standard MATLAB function to add a subtitle to a figure
%    title: Standard MATLAB function to add a title to a figure
%    xlabel: Standard MATLAB function to label the X-axis of a plot
%    ylabel: Standard MATLAB function to label the Y-axis of a plot
%
%  START OF EXECUTABLE CODE
%
%  Create an array of all harmonics included in the input data
%  (including DC and fundamental)
%
   exciteFreq = freq * ( 0:( nHarmonics - 1 ) );
%
%  Single out the highest harmonic included in the input data
%
   cutoffFreq = exciteFreq( nHarmonics );
%
%  Create a new figure window for the spectrum plot
%
   figure( figNo );
%
%  Create the spectrum plot to include
%  DC thru highest available harmonic
%
   stem( exciteFreq, harmMag, 'b' );
%
%  Scale the x- and y-axes
%
   axis( [ 0, cutoffFreq, 0, 1.2 * max( harmMag ) ] );
%
%  Label the x- and y-axes
%
   xlabel( 'Frequency (Hz)', 'interpreter', 'latex' );
   ylabel( 'Magnitude (N)', 'interpreter', 'latex' );
%
%  Add a title
%
   title( 'Frequency Content of Input Pulse', ...
          'interpreter', 'latex', 'fontsize', 16 );
%
%  Convert key numerical parameters of the analysis 
%  to printable text
%
   pulseWidthTxt = num2str( pulseWidth );
   subtitleText = [ 'Pulse Width: ', pulseWidthTxt, ' s' ]; 
   fundamentalTxt = num2str( freq );
   harmonicsTxt = num2str( nHarmonics );
%
%  Label the critical parameters of the analysis in a subtitle
%
   subtitleText = [ 'Pulse Width: ', pulseWidthTxt, ' s', ...
                    '\hspace{2.0em} $F_0$: ', ...
                    fundamentalTxt, ' Hz \hspace{2.0ex}', ...
                    harmonicsTxt, ' Harmonics' ];
   subtitle( subtitleText, 'interpreter', 'latex' );
end








function [int_y] = simpson34( x_min, x_max, numPoints, r, forceArray)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME: simpson34
%
% PURPOSE: To integrate a function using Simpson's 1/3 rule and Simpson's
% 3/8 rule (when required) over a given interval using a given number of
% data points.
%
% INPUT: lower integral bound, upper integral bound, number of data points,
% equation to be integrated
%
% OUTPUT: integral of equation from lower to upper integral bound using
% Simpson's 1/3 and 3/8 rule (when required)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Sophia Steblay
% DATE: 10/10/2024
%
% DESCRIPTION OF LOCAL VARIABLES: 
%
% num13Rule: number of points left to integrate between after using 
% Simpson's 3/8 rule
% 
% x_end38: last x-value used in Simpson's 3/8 rule calculation
%
% int_yStore: temporary variable used to store int_y from Simpson's 3/8
% rule
%
% FUNCTIONS CALLED: None
%
% START OF EXECUTABLE CODE
%

    % Calculate h value (aka step)
    step = (x_max - x_min) / (numPoints - 1);
    x = zeros(length(r));
   % x=linspace(r(1),r(end),length(r));
   % y = zeros(length(r));
    x(1) = r(1);
    y=forceArray;

    
   % for i = 2:1:length(r)
    %    x(i) = step * i;
       % y(i) = interp1(r, forceArray, x(i));
   %  end
    
    % Initialize integral of equation to 0
    int_y = 0;
   
    % If the desired number of data points is even, use Simpson's 3/8 rule
    % once and use Simpson's 1/3 rule for the remaining data points
    if rem(numPoints, 2) == 0

        % Find the integral using the formula for Simpson's 3/8 rule
        int_y = int_y + y(1) + y(4);
        int_y = int_y + 3 * y(2);
        int_y = int_y + 3 * y(3);
        int_y = int_y * step * 3/8;

        % Store this integral result in a temporary variable
        int_yStore = int_y;

        % Begin calculations for Simpson's 1/3 rule from the data point
        % that Simpson's 3/8 rule ended with
        int_y = y(5) + y(end);
      
        % Loop through the remaining data points after Simpson's 3/8 rule
        % was applied. Calculate the integral using Simpson's 1/3 rule
        for i = 6:1:(numPoints - 1)

            % if the data point is in an even position
            if rem(i, 2) == 0
                int_y = int_y + 2 * y(i);
              
            else
                int_y = int_y + 4 * y(i);
              
            end

        end

        int_y = int_y * step / 3;
       
        % Add integral result from Simpson's 3/8 rule to the result from
        % Simpson's 1/3 rule to get final integral result
        int_y = int_y + int_yStore;
       

    % If number of data points is odd, use Simpson's 1/3 rule
    else
        
        % Begin integrating using Simpson's 1/3 rule
        int_y = int_y + y(1) + y(end);

        % Iterate through the desired number of data points
        for i = 2:1:(numPoints)

            % if the data point is in an even position
            if rem(i, 2) == 0
                int_y = int_y + 2 * y(i);
              
            else
                int_y = int_y + 4 * y(i);
              
            end

        end

        % Find final integral result using Simpson's 1/3 rule
        int_y = int_y * step / 3;
    end

end % End function








function [max_CP, optimal_pitch, optimal_lambda, lambda_values, pitch_angles, CP_matrix] = sweep_pitch_and_lambda(U, lambda_min, lambda_max, lambda_step, pitch_min, pitch_max, pitch_step, axialFactorMax, angularFactorMax, r, chord, twist, airfoil_types, airfoil_perf, rho, mu, dr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:
%   sweep_pitch_and_lambda
%
% PURPOSE:
%   Sweep through ranges of blade pitch angles and tip speed ratios to
%   determine the combination that yields the max coefficient of power for
%   a specified wind velocity.
%
% INPUT:
%   U - Wind velocity [m/s]
%   lambda_min - Minimum tip speed ratio
%   lambda_max - Maximum tip speed ratio
%   lambda_step - Step size for tip speed ratio sweep
%   pitch_min - Minimum pitch angle [degree]
%   pitch_max - Maximum pitch angle [degree]
%   pitch_step - Step size for pitch angle sweep [degree]
%   axialFactorMax - Axial induction factor
%   angularFactorMax - Angular induction factor
%   r - Radial positions along the blade [m]
%   chord - Chord lengths at radial positions [m]
%   twist - Blade twist angles at radial positions [degree]
%   airfoil_types - Array of airfoil types corresponding to radial positions
%   airfoil_perf - Structure containing airfoil performance data for various types
%   rho - Air density [kg/m^3]
%   mu - Dynamic viscosity of air [kg/(m.s)]
%   dr - Differential segment lengths along the blade [m]
%
% OUTPUT:
%   max_CP - Maximum coefficient of power achieved
%   optimal_pitch - Blade pitch angle corresponding to max_CP [degree]
%   optimal_lambda - Tip speed ratio corresponding to max_CP
%   CP_matrix - Matrix of CP values [pitch_angles x lambda_values]    
%   pitch_angles - Array of swept pitch angles [degree]                    
%   lambda_values - Array of swept tip speed ratios            



% LOCAL VARIABLES :
% lambda_values-array which sweeps the tip speed ratio of the blades 
% pitch-angles- array which ranges from the min to the max pitch angle 
% current_pitch - current pitch angle [°]
% current_lambda - current tip speed ratio
% current_omega - rotational velocity based on the current tip speed ratio
% CP_Temp- return value of calculate_cp_ct main which is used to find the
% optimal pitch angle and tip speed ratio
% optimal_pitch - max pitch angle which achieved the highest cp value 
% optimal_lambda- max tip speed ratio which achieved the highest cp value 

% Function Called: 
%Calc_Cp_Ct_Main- function which returns the coefficient of power and
%thrust
% max - returns the max value in a array 
% size - finds the size of an array 
% ind2sub - converts linear indices to subscripts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define ranges for pitch angles and tip speed ratios
    lambda_values = lambda_min:lambda_step:lambda_max;
    pitch_angles = pitch_min:pitch_step:pitch_max;

    % Initialize CP matrix
    CP_matrix = zeros(length(pitch_angles), length(lambda_values));

    % Nested loops to sweep through all combinations
    for i = 1:length(pitch_angles)
        for j = 1:length(lambda_values)
            current_pitch = pitch_angles(i);
            current_lambda = lambda_values(j);

            % Calculate rotational velocity (omega) based on tip speed ratio
            current_omega = (current_lambda * U) / 48; % [rad/s], R = 48 m

            % Call calc_CP_CT function
            [CP_temp, ~] = Calc_Cp_Ct_Main(U, axialFactorMax, angularFactorMax, r, chord, airfoil_types, current_omega, current_pitch, rho, mu,twist);
         
            % Check if CP exceedes Betz limit, then store CP in the matrix
            if CP_temp > 0.593 || CP_temp < 0
                CP_matrix(i, j) = 0;
            else
                CP_matrix(i, j) = CP_temp;
            end
        end
    end

    % Find the max CP and its indices
    [max_CP, idx] = max(CP_matrix(:), [], 'omitnan');
    [row, col] = ind2sub(size(CP_matrix), idx);

    % Retrieve the optimal pitch angle and tip speed ratio
    optimal_pitch = pitch_angles(row);
    optimal_lambda = lambda_values(col);

    if isnan(max_CP) || isempty(max_CP)
        warning('No valid CP values found within the Betz limit.');
    end
end













function [U_rel, AoA] = AngleOfAttack_debug(U, alpha, alpha_p, r, omega, theta_P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:
%   AngleOfAttack
%
% PURPOSE:
%   Computes the relative wind speed magnitude (U_rel) and the angle of
%   attack (AoA) for a horizontal-axis wind turbine (HAWT) blade section.
%
% INPUT:
%   U - Wind velocity [m/s]
%   alpha - Axial induction factor [unitless]
%   alpha_p - Angular induction factor [unitless]
%   r - Radial position [m]
%   omega - Rotational velocity [rad/s]
%   theta_P - Section pitch angle [degrees]
%
% OUTPUT:
%   U_rel - Relative velocity magnitude [m/s]
%   AoA - Angle of attack [degrees]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if length(theta_p)==1
%
U_axial = U * (1 - alpha); % Compute the axial component of the wind speed
U_tangential = omega * r * (1 + alpha_p); % Compute the tangential component of the wind speed
U_rel = sqrt(U_axial^2 + U_tangential.^2); % Compute the relative wind speed magnitude
psi = atan2(U_axial, U_tangential); % Compute the angle of the relative wind (phi) in radians
theta_P_rad = deg2rad(theta_P); % Convert the section pitch angle to radians
AoA_rad = psi - theta_P_rad; % Compute the angle of attack in radians
AoA=AoA_rad; % DEBUG - changed line - Convert the angle of attack to degrees

end