%---------------------------------------------------
%  NAME:      example_3D_solenoid_filament.m
%  WHAT:      Calculation of the magnetic field of a solenoid
%             on a volume (+ 3D plot).
%  REQUIRED:  BSmag Toolbox 20150407
%  AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%----------------------------------------------------

% Initialize
clear all, close all, clc
BSmag = BSmag_init(); % Initialize BSmag analysis

% Source points (where there is a current source)
theta = linspace(-2*2*pi,2*2*pi,200);
Gamma = [cos(theta'),sin(theta'),theta'/10]; % x,y,z [m,m,m]
%vorm van het figuur

I = 1; % filament current [A]
dGamma = 1e9; % filament max discretization step [m]
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

% Field points (where we want to calculate the)
x_M = linspace(-1.5,1.5,10); % x [m]
y_M = linspace(-2,2,10); % y [m]
z_M = linspace(-2,2,10); % z [m]
[X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);

%BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points volume

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);

% Plot B/|B|
figure(1)
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')
%axis tight