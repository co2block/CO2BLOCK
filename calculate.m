function [d_list,well_list,d_max,Q_M_each,V_M] = calculate(fpath,fname,correction,dist_min,...
    dist_max,nr_dist,nr_well_max,rw,time_yr,maxQ)

    %read data
    [thick,area_res,perm,por,dens_c,visc_c,visc_w,compr,p_lim,rc,gamma,delta, omega] = read_data(fpath,fname);
    time = time_yr*86400*365 ;                                              %injection time [sec]
    R_influence = sqrt(2.246*perm*time/(visc_w*compr));                     % pressure propagation radius for the time of injection

    % initialize
    if nr_well_max == 'auto'                                                % calculate maximum well number if not set
        nr_well_max = floor(area_res/(dist_min^2));
    end
    
    if dist_max == 'auto'                                                   % calculate maximum interwell distance if not set
        dist_max = sqrt(area_res)/2;
    end
    
    d_list = linspace(dist_min,dist_max,nr_dist) ;                          % inter-well distance list
    M0 = 1 ;                                                                % guess value for total injection rate [Mton/y] - any value


    b = []; 
    well_list = [];   
    d_max = [];   
    w_id = 0 ;

    for x_grid_num = 1:sqrt(nr_well_max)                                    % number of wells on a horizontal row
        for y_grid_num = x_grid_num:x_grid_num+1                            % number of wells on a  vertical row
            w_id = w_id +1 ;                                                % well number scenario ID
            w = x_grid_num*y_grid_num ;                                     % well number for each scenario
            well_list(w_id) = w;                                            % store in vector
            d_max(w_id) = sqrt(area_res/w);                                 % maximum interwell distance for each well number [km]
            for d = 1:nr_dist
                % calculate grid and distances
                distance= d_list(d)*1000;
                wells_coord_x = repmat((0:distance:distance*x_grid_num-1),[y_grid_num,1]) ;
                wells_coord_y = repmat(transpose(0:distance: distance*y_grid_num-1),[1,x_grid_num]) ;
                central_well_x = ceil(x_grid_num/2);                             % posiiton of the central well in x-coord vector
                central_well_y =  ceil(y_grid_num/2);                            % posiiton of the central well in y-coord vector
                dist_vec_x = wells_coord_x - wells_coord_x(central_well_y,central_well_x);        % distance in x from central well [km]
                dist_vec_y = wells_coord_y - wells_coord_y(central_well_y,central_well_x) ;       % distance in y from central well [km]
                dist_vec   = sqrt(dist_vec_x.^2+dist_vec_y.^2) ;                                  % distance from central well [m]
                dist_vec(central_well_y,central_well_x) = rw ;                                    % assign wells radius to the central well

                % calculate pressure build-up for a reference flow rate Q0 at each scenario
                Q0 = M0*1e9/dens_c/365/86400/w;                                % injection rate per well [m3/s]
                Q0_vec(w_id) = Q0;                                             % store in vector
                csi = sqrt(Q0*time/pi/por/thick);                              % average plume extension [m]
                psi = exp(omega)*csi;                                          % equivalent plume extension [m]
                p_c = (Q0*visc_w)/(2*pi*thick*perm)/1e6;                       % characteristic pressure  [MPa]
                d_min_p(w_id) = 2*csi/1000 ;                                   % minimum interwell distance for each well number case [km]
                p_sup = 0; 
                for i = 1:w
                    r = dist_vec(i);
                    Delta_p = Nordbotten_solution(r,R_influence,psi,rc,gamma)*p_c;       % overpressure according to Nordbotten and Celia solution for overpressure [MPa]
                    p_sup = p_sup +  Delta_p ;                                 % superposed overpressure [MPa]
                end

                switch correction                                              % correction for superposition error (De Simone et al., GRL2019)
                    case 'off'
                        prefix = '';
                        sup_error = 0 ;
                        b(w_id) = (visc_w-visc_c)/4/pi/perm/thick  ;                   
                    case 'on'  
                        prefix = 'corrected_' ;
                        if w < 9 ||  R_influence*csi/(distance^2) < 1
                            sup_error = 0; 
                            b(w_id) = (visc_w-visc_c)/4/pi/perm/thick  ;
                        else      
                            sup_error = w*delta/4 * log(R_influence*csi/(distance^2)) ;  
                            b(w_id) = (visc_w-visc_c)/4/pi/perm/thick * (1+w/4) ; 
                        end
                end
                p_sup =  p_sup - sup_error*p_c ; 
                p_sup_vec(w_id,d) =  p_sup  ;          
            end
        end
    end  

    % calculate injectable per well flow-rate Q_M_each and total storage volume V_M for each scenario
    b = repmat(b',1,nr_dist,1);   
    q1 = repmat(Q0_vec',1,nr_dist); 
    well_mat = repmat(well_list',1,nr_dist); 

    p1 = p_sup_vec*1e6;
    p2 = p_lim*1e6;   
    q2 = - p2./b./(lambertw(-1,-p2./q1./b.*exp(-p1./q1./b))) ;     % [m3/s] limit flow rate at each well with  non-linear multi-phase relationship
    Q_M_each = q2*86400*365*dens_c/1e9 ;                    %[ Mton/year]   limit flow rate at each well with non-linear multi-phase relationship 

    %rescale according to lower contstraints
    for dd = 1:nr_dist
        distance= d_list(dd)*1000;
        if Q_M_each(1,dd) > 0.9999*maxQ            % rescale for n_well = 1
           Q_M_each(1,dd) =  0.9999*maxQ ;
        end
        for nn = 2:w_id                             % rescale for n_well > 1
            if Q_M_each(nn,dd) > 0.9999*maxQ || q2(nn,dd) > (distance)^2*pi*por*thick/4.001/time 
               Q_M_each(nn,dd) = min( 0.9999*maxQ, (distance)^2*pi*por*thick/4.001/time*86400*365*dens_c/1e9) ;
            end
        end
    end


    Q_M_tot = Q_M_each.*well_mat ;                     %[ Mton/year]   total sustainable flow rate   
    V_M = Q_M_tot.*time_yr/1000 ;                      %[Gton] total sustainable injectable volume


    % upper constraint 
    d_mat = repmat(d_list,w_id,1); 
    d_max = repmat(d_max,nr_dist,1); 
    d_max_check =  d_mat./d_max' ;       

    % apply upper constraint
    possible = d_max_check < 1 ; 
    V_poss = V_M .* possible;
    well_poss = well_mat.* possible;
    d_poss = d_mat.* possible;

end