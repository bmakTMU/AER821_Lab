% SETUP_SRV02_EXP1
%
% SRV02-ET Experiment - Position Control
% 
% SETUP_SRV02_EXP1 sets the SRV02 system model parameters accordingly to 
% the user-defined configuration.
% SETUP_SRV02_EXP1 also sets the controller parameters, accordingly to 
% the user-defined specifications.
%
% Copyright (C) 2002 Quanser Consulting Inc.
% Quanser Consulting Inc.

clear all;


% ############### USER-DEFINED SRV02 CONFIGURATION #####################################
% External Gear Configuration: set to 'HIGH' or 'LOW'
EXT_GEAR_CONFIG = 'HIGH';
% Encoder Type: set to 'E' or 'EHR'
ENCODER_TYPE = 'E';
% Is the SRV02 equipped with a Tachometer? (i.e. option T): set to 'YES' or 'NO'
TACH_OPTION = 'YES';
% Universal Power Module (UPM) Type: set to 'UPM_2405', 'UPM_1503', or 'UPM_1503x2'
UPM_TYPE = 'UPM_1503';
% Type of Load: set to 'NO_LOAD', 'DISC_LOAD', or 'BAR_LOAD'
LOAD_TYPE = 'DISC_LOAD';
% Cable Gain used: set to 1, 3, or 5
K_CABLE = 1;
% ############### END OF USER-DEFINED SRV02 CONFIGURATION ###############################


% Set Model Variables Accordingly to the USER-DEFINED SRV02 System Configuration
% Also Calculate the SRV02 Model Parameters and 
% Write them to the MATLAB Workspace (to be used in Simulink diagrams)
[ Rm, Kt, Km, Kg, Eff_G, Beq, Jeq, Eff_M ] = Set_SRV02_Configuration ( EXT_GEAR_CONFIG, ENCODER_TYPE, TACH_OPTION, UPM_TYPE,LOAD_TYPE );

% Define which Global Variables will be required in Simulink diagrams
global K_POT K_TACH K_ENC VMAX_UPM IMAX_UPM

% ############### USER-DEFINED CONTROLLER SPECIFICATIONS ###############################
% Controller Parameters
zeta = sqrt(2)/2; % 5% overshoot
Tp = 0.150; % 150 ms time of first peak

% Student Configuration: Set to 'MANUAL' if you want the student to calculate gains
%                        or set to 'AUTO' if you would like this script to set the gains 
%                        to meet the desired specifications
Student_Config = 'MANUAL';
% ############### END OF USER-DEFINED CONTROLLER SPECIFICATIONS ########################


if strcmp ( Student_Config, 'AUTO' )
    % Design the controller as per the specifications above
    [ Kp, Kv ] = d_SRV02_Position_PV ( Rm, Kt, Km, Kg, Eff_G, Beq, Jeq, Eff_M, zeta, Tp );
    % Display the calculated gains
    disp( 'PV controller gains: ' )
    disp( [ 'Kp = ' num2str( Kp ) ' V/rd' ] )
    disp( [ 'Kv = ' num2str( Kv ) ' V/rd/s' ] )
else
    disp( '############################################' )
    disp( 'Kp and Kv (Controller Gains) need to be set.' )
    disp( '############################################' )
    Kp = 0, Kv = 0 
end