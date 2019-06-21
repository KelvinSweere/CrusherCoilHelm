    %---------------------------------------------------
    %  CRUSHER COIL HELM - kopie example_3D_solenoid_filament.m
    %  Kelvin Sweere
    %  Voor het instaleren moet de BSmag Core aan path worden
    %  toegevoegd.
    %----------------------------------------------------

    % Initialize
    clear all, close all,  clc
    %init toolbox
    BSmag = BSmag_init(); % Initialize BSmag analysis

    % Hoeveel resolutie per as 
    x_res = 100;
    y_res = 100;
    z_res = 100;
    %speling vanaf assen
    fct_z = 0.1;
    fct= 0.1;
    %stroom door draad
    I = 10.0; % filament current [A]
    %wikkeling frequentie
    freq = 20;
    %Helm in het figuur weergeven.
    helm_aan = 1;   %aan!
   
    %% vorm van het figuur
    t = 0.5*pi*(3/18):pi/2000:pi/2;         %tijd
    xt(1,:) = 3*sin(0.8*t).*(0.18*cos(freq*t));
    yt(1,:) = 3*sin(0.8*t).*(0.18*sin(freq*t)*1.1164);
    z = 0.6*cos(t);
    
    
    %% tweede figuur
    xt_tmp(1,:) = 3*sin(0.8*t).*(0.18*cos(freq*t+(pi/3)));
    yt_tmp(1,:) = 3*sin(0.8*t).*(0.18*sin(freq*t+(pi/3)*1.1164));
    z_tmp = 0.6*cos(t);
    
    % gamma waardes instellen
    Gamma_1 = [xt',yt',z']; % x,y,z [m,m,m]
    Gamma_2 = [xt_tmp',yt_tmp',z_tmp']; % x,y,z [m,m,m]
    
    % Opdelen in kleinere stukjes
    dGamma = 1e9; % filament max discretization step [m]
    [BSmag] = BSmag_add_filament(BSmag,Gamma_1,I,dGamma);
    [BSmag] = BSmag_add_filament(BSmag,Gamma_2,-1*I,dGamma);
    
    %% code na de lange berekeningen
    %plot het tweede gemaakte figuur ook in hetzelfde figuur.
    plot3(xt,yt,z);
    plot3(xt_tmp,yt_tmp,z_tmp);    
    
    % Maken van een meshgrid x,y,z
    x_M = linspace(min(xt)-fct,max(xt)+fct,x_res); % x [m]
    y_M = linspace(min(yt)-fct,max(yt)+fct,y_res); % y [m]
    z_M = linspace(min(z)-fct_z,max(z)+fct_z,z_res); % z [m]
    [X_M,Y_M,Z_M] = meshgrid(x_M,y_M,z_M);

    %Ik heb deze uitgezet omdat ik dit niet echt duidelijk vindt.
    %BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points volume

    % Biot-Savart Integration
    [BSmag_1,X,Y,Z,BX,BY,BZ]        = BSmag_get_B(BSmag,X_M,Y_M,Z_M);
    %tweede figuur
    [BSmag_2,X,Y,Z,BX_2,BY_2,BZ_2]  = BSmag_get_B(BSmag,X_M,Y_M,Z_M);

%%     Plot B/|B|
    figure(1)
        normB=sqrt(BX.^2+BY.^2+BZ.^2);
%         teken de vectoren
        quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b');
    axis tight

%% Plot Bz on the volume
    figure(2), subplot(2,2,1),hold on, box on, grid on
    if helm_aan == 1 
        plot3(Gamma_1(:,1),Gamma_1(:,2),Gamma_1(:,3),'.-r'); % plot filament 
        plot3(Gamma_2(:,1),Gamma_2(:,2),Gamma_2(:,3),'.-b'); % plot filament    
    end
    %xz
    s = slice(X,Y,Z,BZ_2,[],[0],[]); % plot B
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    s.DiffuseStrength = 0.8;
    colorbar;
    colormap('jet');    %default
    caxis([-5,5]*1e-5);
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]');, title ('Bz [T]');
    view(0,0), axis equal, axis tight, legend   %  caxis([-Tesla,Tesla]*1e-5)      %bepaal max- en minimale waardes voor kleuren

    %%  Figuur 2 (dwarsdoorsnedes)
    subplot(2,2,2), hold on, box on, grid on
    if helm_aan == 1 
        plot3(Gamma_1(:,1),Gamma_1(:,2),Gamma_1(:,3),'.-r'); % plot filament
        plot3(Gamma_2(:,1),Gamma_2(:,2),Gamma_2(:,3),'.-b'); % plot filament
    end
    s = slice(X,Y,Z,BZ,[0],[],[]); % plot Bz
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    %s.DiffuseStrength = 0.8;
    colorbar;
    colormap('jet');    %default
    caxis([-5,5]*1e-5);
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
    view(90,0), axis equal, axis tight, legend   %bepaal max- en minimale waardes voor kleuren

    %     3
    subplot(2,2,3),hold on, box on, grid on
    if helm_aan == 1 
        plot3(Gamma_1(:,1),Gamma_1(:,2),Gamma_1(:,3),'.-r'); % plot filament    
        plot3(Gamma_2(:,1),Gamma_2(:,2),Gamma_2(:,3),'.-b'); % plot filament    
    end
    s = slice(X,Y,Z,BZ,[],[],[0.35]); % plot Bz   %was 0.3
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    colorbar;
    colormap('jet');    %default
    caxis([-5,5]*1e-5);
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
    view(), axis equal, axis tight, legend %bepaal max- en minimale waardes voor kleuren
    
%     %Hieronder kan handig zijn maar heb ik niet verder gebruikt.    

% %%     Plot some flux tubes
%     figure(3), hold on, box on, grid on
%     	plot3(Gamma_1(:,1),Gamma_1(:,2),Gamma_1(:,3),'.-r') % plot filament
%     	[X,Y,Z] = ndgrid(-2:0.5:2,-2:0.5:2,-2); % define tubes starting point        
%     	htubes = streamtube(stream3(X,Y,Z,BX,BY,BZ,X0,Y0,Z0), [0.2 10]);
%     xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Flux tubes')
%     view(3), axis equal, axis tight
%     set(htubes,'EdgeColor','none','FaceColor','c') % change tube color
%     camlight left % change tube light
