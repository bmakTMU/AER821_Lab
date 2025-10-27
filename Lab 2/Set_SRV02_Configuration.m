% SET_SRV02_CONFIGURATION
%
% SET_SRV02_CONFIGURATION accepts the user-defined configuration 
% of the Quanser SRV02 system. SET_SRV02_CONFIGURATION then sets up 
% the SRV02 configuration-dependent model variables accordingly,
% and finally returns the calculated model parameters of the SRV02-ET Quanser plant.
%
% SRV02 system nomenclature (and SET_SRV02_MODEL_PARAMETERS returned values):
% Rm        Motor Armature Resistance               (Ohm)
% Kt        Motor Torque Constant                   (N.m/A)
% Km        Motor Back-EMF Constant                 (V.s/rd)
% Kg        Total Gear Ratio
% Eff_G     Gearbox Efficiency
% Eff_M     Motor Efficiency
% Beq       Equivalent Viscous Damping Coefficient 
%                       as seen at the Load         (N.m.s/rd)
% Jeq       Equivalent Inertia
%                       as seen at the Load         (kg.m^2)
% K_POT     Potentiometer Sensitivity               (rd/V)
% K_TACH    Tachometer Sensitivity                  (rd/s/V)
% K_ENC     Encoder Resolution                      (rd/count)
% VMAX_UPM  UPM Maximum Output Voltage              (V)
% IMAX_UPM  UPM Maximum Output Current              (A)
% K_R2D     Conversion Factor: from Radian to Degree
%
% Copyright (C) 2002 Quanser Consulting Inc.
% Quanser Consulting Inc.


% USER-DEFINED SRV02 System Configuration
function [ Rm, Kt, Km, Kg, Eff_G, Beq, Jeq, Eff_M ] = Set_SRV02_Configuration( EXT_GEAR_CONFIG, ENCODER_TYPE, TACH_OPTION, UPM_TYPE, LOAD_TYPE )
% Calculate Useful Conversion Factors
Calc_Conversion_Constants;
% Calculate SRV02 Model Parameters
[ Rm, Kt, Km, Kg, Eff_G, Beq, Jeq, Eff_M ] = Calc_SRV02_Parameters( EXT_GEAR_CONFIG, ENCODER_TYPE, TACH_OPTION, UPM_TYPE, LOAD_TYPE );
% end of 'Set_SRV02_Configuration( )'


% Calculate the Global Variables used in Simulink Diagrams
function [ Rm, Kt, Km, Kg, Eff_G, Beq, Jeq, Eff_M ] = Calc_SRV02_Parameters( EXT_GEAR_CONFIG, ENCODER_TYPE, TACH_OPTION, UPM_TYPE, LOAD_TYPE )
% Variables Set (These variables are needed in the workspace and so must be Global)
global K_POT K_TACH K_ENC VMAX_UPM IMAX_UPM 
% This functions sets the Constant parameters associated with the SRV02
[ Rm, Kt, Km, Kgi, Eff_G, Eff_M ] = Set_SRV02_Constants;
% This function Calculates the inertias present in the system
[ Jmotor, J24, J72, J120 ] = Calc_SRV02_Inertias( TACH_OPTION );

% External Gear Ratio
if strcmp (EXT_GEAR_CONFIG, 'LOW')
    % Low Gear Configuration: (3x) 72-tooth gears
    Kge = 1;
    Kg = Kgi * Kge;
    Jeq = 3 * J72 + Kg^2 * Jmotor * Eff_G;
    % Equivalent Viscous Damping Coefficient as seen at the Load (N.m.s/rd)
    Beq = 1.5e-3;
elseif strcmp (EXT_GEAR_CONFIG, 'HIGH')
    % High Gear Configuration: (1x) 24-tooth gear, (2x) 72-tooth gear, (1x) 120-tooth gear
    Kge = 5;
    Kg = Kgi * Kge;
    Jeq = J24 + 2 * J72 + J120 + Kg^2 * Jmotor * Eff_G; %without backlash gear
    %Jeq = J24 + 2 * J72 + J120 + Kg^2 * Jmotor * Eff_G;
    % Equivalent Viscous Damping Coefficient as seen at the Load (N.m.s/rd)
    Beq = 4e-3;
else
    disp( 'Error: Please Set the SRV02 Gear Configuration.' )
    Kge = 0;
    Kg = 0;
    Jeq = 0;
    Beq = 0;
end

% Set the SRV02 Load Inertia
if strcmp (LOAD_TYPE, 'NO_LOAD')
    Jload = 0;
elseif strcmp (LOAD_TYPE, 'DISC_LOAD')
    Jload = 1/2 * 0.04 * (0.05)^2; % J = mr^2 / 2
elseif strcmp (LOAD_TYPE, 'BAR_LOAD')
    Jload = 1/12 * 0.038 * (0.1525)^2 + 0.038 * (0.0575)^2; % J = mL^2 / 12 + mR^2 (Moment about center then shifted)
elseif strcmp (LOAD_TYPE, 'SLIP_RING')
    mtracker = 0.244; % mass (kg)
    rtracker = 2 * 0.0254; % radius (m)
    Jload = mtracker * rtracker^2 / 2; % J = mr^2 / 2
elseif strcmp (LOAD_TYPE, 'ROTFLEX')
    Jload = 0.0005; % Calculated Experimentally
elseif strcmp (LOAD_TYPE, 'FLEXGAGE')
    Jload = 0.0002; % Calculated Experimentally
elseif strcmp (LOAD_TYPE, 'ROTPEN')
    Jload = 0.0015; % Calculated Experimentally
else
    disp( 'Error: Please Set the SRV02 Load Configuration.' )
    Jload = 0;
end
% Equivalent Inertia of the SRV02 System with Load (kg.m^2)
Jeq = Jeq + Jload;

% Potentiometer Sensitivity (rd/V)
K_POT = -(352 * pi / 180 / 10);
% Tachometer Sensitivity (rd/s/V)
K_TACH = -(1000 * 2 * pi / 60 / 1.5 / Kg);
% Encoder Resolution, for a quadrature encoder, (rd/count)
if strcmp (ENCODER_TYPE,'E')
    K_ENC = 2 * pi / ( 4 * 1024 );
elseif strcmp (ENCODER_TYPE,'EHR')
    K_ENC = 2 * pi / ( 4 * 2048 );
else
    disp( 'Error: Please Set the SRV02 Encoder Type.' )
    K_ENC = 0;
end

% Set the UPM Maximum Output Voltage (V) and Output Current (A)
if  strcmp (UPM_TYPE, 'UPM_2405')
    VMAX_UPM = 24;
    IMAX_UPM = 5;
elseif ( strcmp (UPM_TYPE, 'UPM_1503') | strcmp (UPM_TYPE, 'UPM_1503x2') ) 
    VMAX_UPM = 15;  
    IMAX_UPM = 3;
else
    disp( 'Error: Please Set the UPM Type.' )
    VMAX_UPM = 0;
    IMAX_UPM = 0;
end
% end of 'Calc_SRV02_Parameters( )'


function [ Rm, Kt, Km, Kgi, Eff_G, Eff_M ] = Set_SRV02_Constants( )
global K_RDPS2RPM K_IN2M K_OZ2N
% Armature Resistance (Ohm)
Rm = 2.6;
% Motor Torque Constant (N.m/A)
Kt = 1.088 * K_OZ2N * K_IN2M; % = .00767
% Motor Back-EMF Constant (V.s/rd)
Km = 0.804e-3 * K_RDPS2RPM; % = .00767
% Internal Gear Ratio (of the Planetary Gearbox)
Kgi = 14;
% Gearbox Efficiency
Eff_G =0.9;
% Motor ElectroMechanical Efficiency
Eff_M =0.69;
% end of function 'Set_SRV02_Constants( )'


% Calculate the SRV02 External Gear Inertias (kg.m^2)
function [ Jmotor, J24, J72, J120 ] = Calc_SRV02_Inertias( TACH_OPTION )
global K_OZ2N K_IN2M
% Rotor Inertia (kg.m^2)
Jm = 5.523e-5 * K_OZ2N * K_IN2M; % = 3.9e-7
% Tachometer Armature Inertia, if any (kg.m^2)
if strcmp ( TACH_OPTION, 'YES')
    Jtach = 1e-5 * K_OZ2N * K_IN2M; % = 7e-8
else
    Jtach = 0;
end
% Motor Equivalent Inertia (kg.m^2)
Jmotor = Jm + Jtach;
% External Gears Inertias (kg.m^2)
% J24: 24-tooth Gear Inertia (on the Motor Shaft)
m24 = .005; % mass (kg)
r24 = 0.5 / 2 * 0.0254; % radius (m)
J24 = m24 * r24^2 / 2;
% J72: 72-tooth Gear Inertia (on the Potentiometer Shaft)
m72 = .030; % mass (kg)
r72 = 1.5 / 2 * 0.0254; % radius (m)
J72 = m72 * r72^2 / 2;
% J120: 120-tooth Gear Inertia (on the Load Shaft)
m120 = .083; % mass (kg)
r120 = 2.5 / 2 * 0.0254; % radius (m)
J120 = m120 * r120^2 / 2;
% end of function 'Calc_SRV02_Inertias( )'


% Calculate Useful Conversion Factors
function Calc_Conversion_Constants ()
global K_R2D K_IN2M K_RDPS2RPM K_RPM2RDPS K_OZ2N 
% from radians to degrees
K_R2D = 180/pi;
% from degrees to radians
K_D2R = 1/K_R2D;
% from Inch to Meter
K_IN2M = 0.0254;
% from Meter to Inch
K_M2IN = 1 / K_IN2M;
% from rd/s to RPM
K_RDPS2RPM = 60 / ( 2 * pi );
% from RPM to rd/s
K_RPM2RDPS = 2 * pi / 60;
% from oz-force to N
K_OZ2N = 0.2780139;
% from N to oz-force
K_N2OZ = 1 / K_OZ2N;
% end of 'Calc_Conversion_Constants( )'