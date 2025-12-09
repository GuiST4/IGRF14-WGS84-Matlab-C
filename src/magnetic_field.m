%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Geomagnetic model of the Earth %%%%%%%%%%%%%%%%%%%%
% Inputs = Position (LLA); Year; Gauss Coefficients; Degree of the SH %
%%%% Outputs = B - Magnetic Field Vector (nT); F - Magnitude of B %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B, F] = magnetic_field(height_m, latitude_deg, longitude_deg, year, C_nm, N)

    a = 6371200; % Radius of the Earth in meters
    lon = deg2rad(longitude_deg);

    % Geodetic -> Geocentric
    A = 6378137; % m
    f = 1/298.257223563;
    e = sqrt(f*(2-f));
    Rc = A/sqrt(1-e^2*(sin(deg2rad(latitude_deg)))^2);
    p = (Rc + height_m)*cos(deg2rad(latitude_deg));
    z = (Rc*(1-e^2) + height_m)*sin(deg2rad(latitude_deg));
    r = sqrt(p^2+z^2);
    lat_c = asin(z/r); 
    theta = pi/2 - lat_c;
    epsilon = deg2rad(latitude_deg) - lat_c;

    % Computing Snm & Knm
    [Snm, Knm] = compute_Snm_Knm(N);

    % Computing Pnm & dPnm
    [Pnm, dPnm] = compute_Pnm_dPnm(theta, Knm, N);

    % Gauss Coefficients and Secular Variations
    sz = size(C_nm, 1); 
    g_nm_norm = C_nm(1:sz, 1:sz);
    h_nm_norm = C_nm(1:sz, sz+1:2*sz);       
    g_nm_dot  = C_nm(1:sz, 2*sz+1:3*sz);     
    h_nm_dot  = C_nm(1:sz, 3*sz+1:4*sz);   
    year_table = C_nm(1, end);

    % Time Coefficients
    [g_nm_norm_t, h_nm_norm_t] = time_coeff(g_nm_norm, h_nm_norm, g_nm_dot, h_nm_dot, year_table, year);

    % Unnormalizing the Gauss Coefficients
    [g_nm, h_nm] = unnormalized_coefficients(g_nm_norm_t, h_nm_norm_t, Snm);

    % Computing B
    sumR = 0;
    sumT = 0;
    sumL = 0;

    for n = 1:N
        scale = ((a/r)^(n+2));
        for m = 0:n
            sumR = sumR + scale*(n+1)*Pnm(n+1,m+1)*(g_nm(n+1,m+1)*cos(m*lon) + h_nm(n+1,m+1)*sin(m*lon));
            sumT = sumT + scale*dPnm(n+1,m+1)*(g_nm(n+1,m+1)*cos(m*lon) + h_nm(n+1,m+1)*sin(m*lon));
            sumL = sumL + scale*m*Pnm(n+1,m+1)*(-g_nm(n+1,m+1)*sin(m*lon) + h_nm(n+1,m+1)*cos(m*lon));
        end
    end

    Br = sumR;
    Bt = -sumT;
    Bl = -sumL/sin(theta);
    
    % Geocentric -> NED
    BX = -Bt*cos(epsilon) - Br*sin(epsilon);
    BY = Bl;
    BZ = Bt*sin(epsilon) - Br*cos(epsilon);

    B = [squeeze(BX); squeeze(BY); squeeze(BZ)];
    F = norm(B);

end

function [Snm, Knm] = compute_Snm_Knm(N)

    Snm = zeros(N+1,N+1);
    Knm = zeros(N+1,N+1);
    
    Snm(1,1) = 1;
    Knm(1,1) = 0;
    
    for n = 1:N
        for m = 0:n
            if n == 1
            Knm(n+1,m+1) = 0;
            else 
                Knm(n+1,m+1) = ((n-1)^2-m^2)/((2*n-1)*(2*n-3));
            end
            if m == 0
                Snm(n+1,m+1) = Snm(n,m+1)*(2*n-1)/n;
            else
                if m == 1
                    J = 2;
                else
                    J = 1;
                end
                Snm(n+1,m+1) = Snm(n+1,m)*sqrt((n-m+1)*J/(n+m));
            end
        end
    end

end

function [g_nm, h_nm] = unnormalized_coefficients(g_nm_un, h_nm_un, Snm)

    [m, n] = size(Snm);

    g_nm_un = g_nm_un(1:n,1:m);
    h_nm_un = h_nm_un(1:n,1:m);

    g_nm = Snm.*g_nm_un;
    h_nm = Snm.*h_nm_un;

end

function [Pnm, dPnm] = compute_Pnm_dPnm(theta, Knm, N)

    Pnm = zeros(N+1,N+1);
    Pnm(1,1) = 1;
    
    for n = 1:N
        for m = 0:n
            if m == n
                Pnm(n+1,m+1) = sin(theta)*Pnm(n,m);
            elseif m == n - 1
                    Pnm(n+1,m+1) = cos(theta)*Pnm(n,m+1);
            else
                if m > n - 2
                    Pnm(n-1,m+1) = 0;
                else
                Pnm(n+1,m+1) = cos(theta)*Pnm(n,m+1) - Knm(n+1,m+1)*Pnm(n-1,m+1);
                end
            end
        end
    end
    
    dPnm = zeros(N+1,N+1);
    dPnm(1,1) = 0;
    
    for n = 1:N
        for m = 0:n
            if m == n
                dPnm(n+1,m+1) = sin(theta)*dPnm(n,m) + cos(theta)*Pnm(n,m);
            elseif n == 1 && m == 0
                dPnm(n+1,m+1) = cos(theta)*dPnm(n,m+1) - sin(theta)*Pnm(n,m+1);
            else
                dPnm(n+1,m+1) = cos(theta)*dPnm(n,m+1) - sin(theta)*Pnm(n,m+1) - Knm(n+1,m+1)*dPnm(n-1,m+1);
            end
        end
    end

end

function [g_nm_t, h_nm_t] = time_coeff(g_nm0, h_nm0, g_nm_dot, h_nm_dot, year_table, year)
    
    time = year - year_table;
    g_nm_t = g_nm0 + time.*g_nm_dot;
    h_nm_t = h_nm0 + time.*h_nm_dot;

end