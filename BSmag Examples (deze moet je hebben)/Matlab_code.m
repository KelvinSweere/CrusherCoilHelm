    %---------------------------------------------------
    %  CRUSHER COIL - kopie example_3D_solenoid_filament.m
    %  Auteur: Kelvin Sweere
    %  Note: voor het installeren moet de BSmag Core aan het path worden toegevoegd!
    %  bron: https://nl.mathworks.com/matlabcentral/fileexchange/50434-biot-savart-magnetic-toolbox
    %----------------------------------------------------

    % Initialize
    clear all, close all,  clc
    %init toolbox
    BSmag = BSmag_init(); % Initialize BSmag analysis

    %%init van alle gekozen variabelen
    % Hoeveel resolutie per as
    x_res = 10;
    y_res = 10;
    z_res = 10;
    %extra space vanaf assen
    fct_z = 0.08;
    fct= 0.08;
    %stroom door draad
    I = 1.5; % filament current [A]
    %wikkeling frequentie
    freq = 10;
    afst = 0.025;

    %Constantes
    helm_aan = 1;   %aan!
    aantal_draden = 4;
    Tesla = 5; %-50 en 50 microTesla

    %Maken van een list voor alle x,y,z gegevens
    xt = cell(1,aantal_draden); %maak 4 cellen aan.
    yt = cell(1,aantal_draden); %   ""
    z = cell(1,aantal_draden);

    %% vorm van het figuur
    % = 0.5*pi*(1/18):pi/2000:pi/2;         %tijd
    t = linspace(0,0.2,1000);
    t_1 = t;                                %voor andere figuur.

    afst_tot_nul_z = 0.025/2;

    xt{1} = linspace(0,0.3,1000);
    yt{1} = 0.1*sin(250*t);
    z{1}  = afst_tot_nul_z+0*t;


    %% tweede figuur
    xt{2} = linspace(0,0.3,1000);
    yt{2} = 0.1*sin(250*t+pi/4);
    z{2}  = afst_tot_nul_z+0*t;

    %% derde figuur
    xt{3} = linspace(0,0.3,1000);
    yt{3} = 0.1*sin(250*t);
    z{3}  = -afst_tot_nul_z+0*t;

    %% vierde figuur
    xt{4} = linspace(0,0.3,1000);
    yt{4} = 0.1*sin(250*t+pi/4);
    z{4}  = -afst_tot_nul_z+0*t;

    %% gamma waardes instellen
    Gamma = cell(1,aantal_draden);
    for i = 1:aantal_draden
        Gamma{i} = [xt{i}',yt{i}',z{i}']; % x,y,z [m,m,m]
    end
    % Opdelen in kleinere stukjes
    dGamma = 1e9; % filament max discretization step [m]
    for i = 1:aantal_draden
        if(mod(i,2) == 0)
            I = I * -1;
        end
        [BSmag] = BSmag_add_filament(BSmag,Gamma{i},I,dGamma);
    end

    %% code na de lange berekeningen
    %plot het tweede gemaakte figuur ook in hetzelfde figuur.
    for i = 1:aantal_draden
        plot3(xt{i},yt{i},z{i});
    end

    % Maken van een meshgrid x,y,z
    x_M = linspace(-fct,max(xt{1})+fct,x_res); % x [m]
    y_M = linspace(-afst-fct,afst+fct,y_res); % y [m]
    z_M = linspace(-fct_z-afst,fct_z+afst,z_res); % z [m]
    [X_M,Y_M,Z_M] = meshgrid(x_M,y_M,z_M);

    % Biot-Savart Integration
    BX = cell(1,aantal_draden);
    BY = cell(1,aantal_draden);
    BZ = cell(1,aantal_draden);

    for i = 1:aantal_draden
        [BSmag,X,Y,Z,BX{i},BY{i},BZ{i}]  = BSmag_get_B(BSmag,X_M,Y_M,Z_M);
    end
%%     Plot B/|B|
    figure(1)
        for i = 1:aantal_draden
            normB=sqrt(BX{i}.^2+BY{i}.^2+BZ{i}.^2);
            %   teken de vectoren
            quiver3(X,Y,Z,BX{i}./normB,BY{i}./normB,BZ{i}./normB,'b');
        end
    axis tight

%% Plot Bz on the volume
    figure(2), subplot(2,2,1),hold on, box on, grid on
    if helm_aan == 1
        for i = 1:aantal_draden
            plot3(Gamma{i}(:,1),Gamma{i}(:,2),Gamma{i}(:,3),'.-r'); % plot filament
        end
    end

    %xz
    b_tot = 0;
    for i = 1:aantal_draden
        b_tot= b_tot + BX{i};
    end

    s = slice(X,Y,Z,b_tot,[],[0],[]); % plot B
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    s.DiffuseStrength = 0.8;
    colorbar;
    colormap('jet');    %default
    caxis([-Tesla,Tesla]*1e-5);
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]');, title ('Bz [T]');
    view(0,0), axis equal, axis tight, legend   %  caxis([-Tesla,Tesla]*1e-5)      %bepaal max- en minimale waardes voor kleuren

    %%  Figuur 2 (dwarsdoorsnedes)
    subplot(2,2,2), hold on, box on, grid on
    if helm_aan == 1
        for i = 1:aantal_draden
            plot3(Gamma{i}(:,1),Gamma{i}(:,2),Gamma{i}(:,3),'.-r'); % plot filament
        end
    end

    b_tot = 0;
    for i = 1:aantal_draden
        b_tot= b_tot + BY{i};
    end
    %0.01/2 is gekozen uit de grafiek 2
    s = slice(X,Y,Z,b_tot,[0.01/2],[],[]); % plot Bz
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    %s.DiffuseStrength = 0.8;
    colorbar;
    colormap('jet');    %default
    caxis([-Tesla,Tesla]*1e-5);
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
    view(90,0), axis equal, axis tight, legend   %bepaal max- en minimale waardes voor kleuren

    %     3
    subplot(2,2,3),hold on, box on, grid on
    if helm_aan == 1
        for i = 1:aantal_draden
            plot3(Gamma{i}(:,1),Gamma{i}(:,2),Gamma{i}(:,3),'.-r'); % plot filament
        end
    end
    b_tot = 0;
    for i = 1:aantal_draden
        b_tot= b_tot + BZ{i};
    end
    s = slice(X,Y,Z,b_tot,[],[],[0]); % plot Bz   %was 0.3
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    colorbar;
    colormap('jet');    %default
    caxis([-Tesla,Tesla]*1e-5);
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
    view(), axis equal, axis tight, legend %bepaal max- en minimale waardes voor kleuren
