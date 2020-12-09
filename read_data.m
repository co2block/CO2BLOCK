% read paramaters and evaluate maximum sustainable pressure
function [thick,area_res,perm,por,dens_c,visc_c,visc_w,compr,p_lim,rc,...
    gamma,delta, omega] = read_data(path,name)
 
    % -- default parameters (used in case they are not provided)
    litho_grad = 23 ;                   % lithostatic gradient [MPa/km]
    hydro_grad = 10 ;                   % hydrostatic gradient [MPa/km]
    temp_grad = 33 ;                    % temperaturegradient [C/km]
    def_k0 = 0.7 ;                      % default stress ratio s3/s1  
    def_friction_angle = 30 ;           % default rock friction angle
    def_cohesion = 0;                   % default rock cohesion  [MPa]
    def_cr = 5*10^-4 ;                  % default rock compressibility [MPa^-1]
    def_cw = 3*10^-4 ;                  % default water compressibility [MPa^-1]
    def_salinity = 180000 ;             % default salinity [ppm]
    % --

    % -- read data
    fullFileName = fullfile(path,name);
    data = readtable(fullFileName); 
    
    domain_type = char(data{1,2});                            % domain confinement
    depth = double(data{1,3});                                % shallowest depth of reservoir [m]
    depth_mean = double(data{1,4});                           % mean depth of reservoir [m]
    thick = double(data{1,5});                                % thickness of reservoir [m]
    area_res = double(data{1,6});                             % area of reservoir [km^2]
    perm = double(data{1,7})*10^-15;                          % intrinsic permeability [m^2]
    por = double(data{1,8});                                  % porosity [-]
    cr = double(data{1,9})/1e6;                               % rock compressibility [1/Pa]
    cw = double(data{1,10})/1e6;                              % water compressibility [1/Pa]
    dens_c = double(data{1,11})*1e3;                          % Density of CO2 [kg/m^3]
    visc_c = double(data{1,12})/1e3;                          % Viscosity of CO2 [Pa.s] 
    visc_w = double(data{1,13})/1e3;                          % Viscosity of water[Pa.s] 
    pres0 = double(data{1,14});                               % pressure at the top of the reservoir [MPa]
    pres0_mean = double(data{1,15});                          % pressure at the centre of the reservoir [MPa]
    T0_mean = double(data{1,16});                             % temperature at the centre of the reservoir [C]
    salinity = double(data{1,17})/1e6;                        % aquifer salinity [ppm/1e6]
    s1_tot = double(data{1,18});                              % total maximum principal stress at the reservoir [MPa]
    stress_ratio = double(data{1,19});                        % ratio of principal effective stresses (s3/s1) [-]
    friction = double(data{1,20});                            % rock friction angle  [deg]
    cohesion = double(data{1,21});                            % rock cohesion coefficient [MPa]
    tens_strength = double(data{1,22});                       % rock tensile strength [MPa]
    % --

    switch domain_type
        case 'Open'
            domain_type = 'open';
            rc = inf;
        case 'open'
            domain_type = 'open';
            rc = inf;
        case 'Closed'
            domain_type = 'closed';
            rc = sqrt(area_res*10^6/pi);    
        case 'closed'
            domain_type = 'closed';
            rc = sqrt(area_res*10^6/pi);
    end  

    %%% calculate some parameters if not given
    if pres0 == 0 ||  isnan(pres0)
        pres0 = hydro_grad*depth/1000 ;
    end

    if pres0_mean == 0 ||  isnan(pres0_mean)
        pres0_mean = hydro_grad*depth_mean/1000 ;
    end

    if T0_mean == 0 ||  isnan(T0_mean)
        T0_mean = temp_grad*depth_mean/1000 + 15 ;
    end

    if s1_tot == 0  ||  isnan(s1_tot)
        s1_tot = litho_grad*depth/1000;
    end
    
    if stress_ratio == 0 ||  isnan(stress_ratio)
        stress_ratio = def_k0 ;
    end

    if friction == 0 ||  isnan(friction)
        friction = def_friction_angle ;
    end
    
    if cohesion == 0 ||  isnan(cohesion)
        cohesion = def_cohesion ;
    end

    if tens_strength== 0 ||  isnan(tens_strength)
        tens_strength = cohesion/2;
    end
                
    if cr == 0 ||  isnan(cr)
        cr = def_cr/1e6;
    end

    if cw == 0 ||  isnan(cw)
        cw = def_cw/1e6;
    end

    if salinity == 0 || isnan(salinity)
        salinity = def_salinity/1e6 ;
    end
    
    if dens_c == 0 ||  isnan(dens_c)
        [~,dens_c,~] = eos(T0_mean, pres0_mean,salinity);
    end

    if visc_c == 0 || isnan(visc_c)
        [~,~,visc_c] = eos(T0_mean, pres0_mean,salinity);
    end

    if visc_w == 0 || isnan(visc_w) 
        [visc_w,~,~] = eos(T0_mean, pres0_mean, salinity) ; 
    end


    %%% calculate some useful parameters 
    s1 = s1_tot - pres0;                                                    % effective maximum principal stress [MPa]
    s3 = stress_ratio*s1 ;                                                  % effective minimum principal stress [MPa]
    theta = (1-sin(deg2rad(friction)))/(1+sin(deg2rad(friction))) ;
    p_lim_shear = (stress_ratio-theta)/(1-theta)*s1 ...
        + cohesion*cos(deg2rad(friction))/sin(deg2rad(friction)) ;          % limit overpressure for shear failure [MPa]
    p_lim_tensile = s3 + tens_strength ;                                    % limit overpressure for tensile failure [MPa]
    p_lim = min(p_lim_shear, p_lim_tensile);                                % limit overpressure  [MPa]

    gamma = visc_c/(visc_w);                                                %Non-dimensional viscosity ratio[-]
    delta = (visc_w-visc_c)/visc_w;  
    omega = (visc_c+visc_w)/(visc_c-visc_w)*log(sqrt(visc_c/visc_w))-1 ;
    compr = cr+por*cw ;                                                     %total compressibility  [1/Pa]
end
