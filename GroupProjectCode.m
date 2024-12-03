vf = 0.45; %volume fraction
vm = 0.55; %matrix fraction

%fiber properties - high strength carbon fiber
fiber_density = 1.8*10^3; %kg per cubic meter
fiber_price = 25.1; %dollars per kg
fiber_E = 225*10^9; %Pa
fiber_poisson = 0.2;
fiber_G = 15*10^9; %Pa

%matrix properties - unfilled epoxy
matrix_density = 1.11 * 10^3; %kg per cubic meter
matrix_price = 2.96; %dollars per kg
matrix_E = 3 * 10^9; %Pa
matrix_poisson = 0.4;
matrix_G = 2.1*10^9; %Pa

%givens and loading conditions
Mx = 1.2*10^3 * 1.1; %N m, D11, load * 1.1 for safety reasons
Ny = 1170 * 10^3 * 1.1; %N, A22, load * 1.1 for safety reasons
l = 10; %m
w = 0.5; %m
bending = 0.05; %strain
tensile = 0.002; %strain

%initialize variables
E1 = vf*fiber_E + vm*matrix_E; %Pa
E2 = vf/fiber_E + vm/matrix_E;
E2 = 1/E2; %Pa
G12 = fiber_G*matrix_G/(vm*fiber_G + vf*matrix_G); %Pa
nu12 = vf*fiber_poisson + vm*matrix_poisson;
nu21 = nu12 * E2 / E1;

%layer identifiers
t1 = 0.00001; %0 degree thickness, m
t2 = 0.00001; %90 degree thickness, m

%required stiffnesses
D11 = Mx/bending; %required stiffness (based on bending and deflection)
A22 = Ny/tensile; %required stiffness (based on tensile stress and deformation)

%Get current A22 and D11
    % Q matrix (material coordinates)
    denom = 1 - nu12 * nu21;
    Q11 = E1 / denom;
    Q12 = nu12 * E2 / denom;
    Q22 = E2 / denom;
    Q66 = G12;
    
    Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
    
    % Ply definitions
    theta = [0 90 0]; % Ply angles in degrees
    thickness = [t1 t2 t1]; % Ply thicknesses
    z = [-sum(thickness)/2, -sum(thickness)/2 + t1, sum(thickness)/2 - t1, sum(thickness)/2]; % Ply boundaries
    
    % Initialize matrices
    A = zeros(3, 3);
    D = zeros(3, 3);
    
    % Calculate A and D matrices
    for i = 1:length(theta)
        % Compute Qbar for the ply
        theta_rad = theta(i) * pi / 180; % Convert angle to radians
        m = cos(theta_rad);
        n = sin(theta_rad);
        T = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, m^2 - n^2];
        Qbar = inv(T) * Q * inv(T)';
        
        % Ply thickness and mid-plane distance
        h_ply = thickness(i);
        z_mid = (z(i) + z(i+1)) / 2;
        
        % Add contributions to A and D matrices
        A = A + Qbar * h_ply;
        D = D + Qbar * (h_ply^3 / 12 + h_ply * z_mid^2);
    end
    
    Anew = A(2,2);
    Dnew = D(1,1);

%optimization parameters
R1 = Dnew/D11; %bending
R2 = Anew/A22; %tensile

while R1 < 1 | R2 < 1;
    % Q matrix (material coordinates)
    denom = 1 - nu12 * nu21;
    Q11 = E1 / denom;
    Q12 = nu12 * E2 / denom;
    Q22 = E2 / denom;
    Q66 = G12;
    
    Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
    
    % Ply definitions
    theta = [0 90 0]; % Ply angles in degrees
    thickness = [t1 t2 t1]; % Ply thicknesses
    z = [-sum(thickness)/2, -sum(thickness)/2 + t1, sum(thickness)/2 - t1, sum(thickness)/2]; % Ply boundaries
    
    % Initialize matrices
    A = zeros(3, 3);
    D = zeros(3, 3);
    
    % Calculate A and D matrices
    for i = 1:length(theta)
        % Compute Qbar for the ply
        theta_rad = theta(i) * pi / 180; % Convert angle to radians
        m = cos(theta_rad);
        n = sin(theta_rad);
        T = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, m^2 - n^2];
        Qbar = inv(T) * Q * inv(T)';
        
        % Ply thickness and mid-plane distance
        h_ply = thickness(i);
        z_mid = (z(i) + z(i+1)) / 2;
        
        % Add contributions to A and D matrices
        A = A + Qbar * h_ply;
        D = D + Qbar * (h_ply^3 / 12 + h_ply * z_mid^2);
    end
    
    Anew = A(2,2);
    Dnew = D(1,1);
    R1 = Dnew/D11; %bending
    R2 = Anew/A22; %tensile

    if R1 < R2;
        t1 = t1 + 0.00001;
    else; 
        t2 = t2 + 0.00001;
    end
end

disp(t1)
disp(t2)

%mass and cost calculations
V = (2*t1 + t2)*l*w;
mass = V*vf*fiber_density + V*vm*matrix_density;
cost = V*vf*fiber_density*fiber_price + V*vm*matrix_density*matrix_price;

disp(mass)
disp(cost)
