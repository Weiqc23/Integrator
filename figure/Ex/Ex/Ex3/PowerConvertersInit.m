% Parameters file for PowerConverters example

%   Pierre Giroux 
%   Copyright 2015-2022 Hydro-Quebec, and The MathWorks, Inc.

%% Sample times
% Ts_PWM = 2.5e-6;           % PWM generators sample time default value(s)
% Ts_Control = 50e-6;        % Control systems sample time default value(s)
% Ts_Power = 2.5e-6;         % Simscape Electrical Specialized Power Systems sample time default value (s)

Ts_PWM = 1e-7;   %2.5e-6    %2.5e-8    % PWM generators sample time default value(s)
Ts_Control = 50e-6;        % Control systems sample time default value(s)
Ts_Power = 1e-7;   %2.5e-6   %2.5e-8  % Simscape Electrical Specialized Power Systems sample time default value (s)

%% Grid parameters
Fnom = 60;               % Nominal system frequency (Hz)
Vnom_grid = 25e3;        % Utility nominal voltage (L-L rms)
Psc_grid = 100e6;        % Short-circuit level (VA)
P_Ld1 = Psc_grid/100;    % Ld1 (W)

%% 1-MVA 3-Level Active Rectifier (main DC supply) 

% DC link:
Pnom_dc_3L = 1e6;                            % Nominal DC link Power (VA)
Vnom_dc_3L = 1000;                           % Nominal DC link voltage (V)
H_3L = 1/Fnom*2;                             % DC link stored energy constant(s) = 2 cycles
Clink_3L = Pnom_dc_3L*H_3L*2 /Vnom_dc_3L^2;  % PM DC link capacitor (F)
Vc_Initial_3L = Vnom_dc_3L/2;                % PM capacitor initial voltage (V)

% Tr1  Transformer:
Pnom_3L = Pnom_dc_3L;         % Transformer nominal power (VA)
Vnom_prim_3L = Vnom_grid;     % Nominal primary voltage (V)
m_nom_3L = 0.8;               % Nominal modulation index for 3-Level rectifier
Vnom_sec_3L = 0.5*Vnom_dc_3L/sqrt(2)*sqrt(3)*m_nom_3L;  % Nominal secondary voltage (V)
Lxfo_3L = 0.10;               % Total Leakage inductance (pu)
Rxfo_3L = 0.10/30;            % Total winding resistance (pu)
Rm_3L = 200;                  % Magnetization resistance (pu)
Lm_3L = 200;                  % Magnetization inductance (pu)

% Filter 1:
Qnom_Filter1=0.15*Pnom_dc_3L; % Nominal reactive power (VA)
Fn_Filter1=33*Fnom;           % Tuning frequency (Hz)
Q_Filter1=5;                  % Quality factor 

%  Control Parameters
Fc_3L=33*Fnom;                % PWM carrier frequency (Hz)
Freq_Filter=1500;             % Measurement filters natural frequency (Hz)

% VDC regulator (VDCreg)
Kp_VDCreg_3L= 2.5;            % Proportional gain
Ki_VDCreg_3L= 275;            % Integral gain
LimitU_VDCreg_3L= 1.5;        % Output (Idref) Upper limit (pu)
LimitL_VDCreg_3L= -1.5;       % Output (Idref) Lower limit (pu)

% Current regulator (Ireg)
Rff_3L= Rxfo_3L;              % Feedforward R
Lff_3L= Lxfo_3L;              % Feedforward L
Kp_Ireg_3L= 0.2/2;            % Proportional gain
Ki_Ireg_3L= 8;                % Integral gain
LimitU_Ireg_3L= 1.5;          % Output (Vdq_conv) Upper limit (pu)
LimitL_Ireg_3L= -1.5;         % Output (Vdq_conv) Lower limit (pu)

%% 2-MVA Twin 2-Level Statcom 

Pnom_stat=2e6;                % Nominal power (VA)                        
Vnom_prim_stat=Vnom_grid;     % Primary nominal voltage (LLrms)
Vnom_sec_stat=1800*1.25;      % Secondary nominal voltage (LLrms)
Lxfo_Tr_stat= 0.06;           % Total Leakage inductance (pu)
Rxfo_Tr_stat= 0.06/30;        % Total winding resistance (pu)
Rm_Tr_stat=200;               % Magnetization resistance (pu)
Lm_Tr_stat=200;               % Magnetization inductance (pu)
Lbase_sec=Vnom_sec_stat^2/Pnom_stat/(2*pi*Fnom);
Lstat=0.1*Lbase_sec;          % phase reactor inductance (H) = 0.1 pu
Rstat=Lstat*(2*pi*Fnom)/100;  % phase reactor resistance (Ohms)

% DC Link
Vnom_dc_stat=2400;            % Nominal DC link voltage (V)
H_stat=1/Fnom*1;              %  DC link stored energy constant(s) = 1 cycles
Clink_stat= Pnom_stat*H_stat*2 /Vnom_dc_stat^2;  %  DC link capacitor (F)
Vc_Initial_stat=Vnom_dc_stat; % Initial DC link capacitor voltage (V)

% Filter
Cfilter=200e-6;
Rfilter=1/(Cfilter*377)/50;

% Control Parameters

Fc_stat=27*Fnom;              % PWM Carrier Frequency (Hz)
Freq_Filter_stat=1500;        % Measurement filter natural frequency (Hz)
RateLimit_Qstat=Fnom*2;       % Limit rising/falling rate of Qref (pu/sec)

% 2-Level rectifier VDC regulator (VDCreg)
Kp_VDCreg_stat=10;            % Proportional gain
Ki_VDCreg_stat= 800;          % Integral gain
LimitU_VDCreg_stat=1.5;       % Output (Idref) Upper limit (pu)
LimitL_VDCreg_stat= -1.5;     % Output (Idref) Lower limit (pu)

% Current regulator (Ireg)
Rff_stat= Rxfo_Tr_stat+0.06/100*2; % Feedforward R
Lff_stat= Lxfo_Tr_stat+0.06*2;     % Feedforward L
Kp_Ireg_stat= 0.35;        % Proportional gain
Ki_Ireg_stat= 80;          % Integral gain
LimitU_Ireg_stat= 1.5;     % Output (Vdq_conv) Upper limit (pu)
LimitL_Ireg_stat= -1.5;    % Output (Vdq_conv) Lower limit (pu)

%% DC Motor Drive

Ra_La=[0.0597,0.0009];     % Armature resistance Ra(Ohms) and inductance La(H)
Rf_Lf=[200,160];           % Field resistance Rf(Ohms) and inductance Lf(H)
Laf=2.621;                 % Field armature mutual inductance (H)
J=10;                      % Total inertia J (kg.m^2)
Bm=0.272;                  % Viscous friction coefficient (N.m.s)
w_Initial=1200*pi/30;      % Initial speed (rad/s)
If_Initial=2.5;            % Initial field current (A)
%
wref_Initial=w_Initial*30/pi; % Initial motor speed (rpm)
Tload=1000;                   % Load torque (N.m)
Fc_motor=2000;                % PWM carrier frequency (Hz)

% Limit rising/falling rate speed setpoint(wref)
RateLimit_wref_motor=2000;    % rpm/s

% Limit rising/falling rate current reference(Iref)
RateLimit_Iref_motor=20000;   % A/s

% Speed regulator (wreg)
Kp_wreg_motor=15;                % Proportional gain
Ki_wreg_motor= 50;               % Integral gain
Limit_wreg_motor= [1000, -1000]; % Output limit [Upper Lower] (Iref)

% Current regulator (Ireg_motor)
Kp_Ireg_motor=0.0015;       % Proportional gain
Ki_Ireg_motor= 0.25;        % Integral gain
Limit_Ireg_motor= [1, 0];   % Output limit [Upper Lower] (D) 

%% 60-Hz Load

Lchoke_HB=500e-6;           % Choke inductance (H)
Rchoke_HB=5e-3;             % Choke resistance (Ohms)
Vnom_Load_HB=340;           % Nominal load voltage (Vrms)
P_Load_HB=400e3;            % Active power (W)
Q_Load_HB= 50e3;            % Capacitive reactive power (var)

Fc_HB=33*Fnom;              % PWM carrier frequency (Hz)
m_HB=1;                     % Modulation index (Uref magnitude)

%% Variable DC Load

Lchoke_Buck=5e-3;           % Choke inductance (H)
Rchoke_Buck=10e-3;          % Choke resistance (Ohms)
Rload_Buck=0.1;             % Load resistance (Ohms)

Vload_Buck=500;             % DC source voltage (V)
AmplitudeVar_Buck=50;       % Vload variation amplitude(V)
FrequencyVar_Buck=5;        % Vload variation frequency (Hz)

Fc_Buck=2000;               % PWM carrier frequency (Hz)
D_Buck=0.55;                % Duty cycle

%%DC Supply

Lchoke_Boost=1e-3;          % Choke inductance (H)
Rchoke_Boost=50e-3;         % Choke resistance (Ohms)
V_DCsrc=500;                % DC source voltage (V)
R_DCsrc=0.5;                % DC resistance (Ohms)
%
Fc_Boost=2000;              % PWM carrier frequency (Hz)
D_Boost=0.8;                % Duty cycle

%% 50-Hz Load

Lchoke_FB=1e-3;             % Choke inductance (H)
Rchoke_FB=50e-3;            % Choke resistance (Ohms)
Fsys_FB=50;                 % Nominal frequency (Hz)
Vnom_Load_FB=600;           % Nominal load voltage (Vrms)
P_Load_FB=200e3;            % Active power (W)
Q_Load_FB= 50e3;            % Capacitive reactive power (var)

Fc_FB=33*Fsys_FB;           % PWM carrier frequency (Hz)
m_FB=0.9;                   % Modulation index (Uref magnitude)
