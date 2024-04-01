clear all;
L = 3;
k = @(t) ((1/L)* ...
    [cosd(t)^2        cosd(t)*sind(t)  -cosd(t)^2       -cosd(t)*sind(t); ...
     cosd(t)*sind(t)  sind(t)^2        -cosd(t)*sind(t) -sind(t)^2      ; ...
     -cosd(t)^2       -cosd(t)*sind(t) cosd(t)^2        cosd(t)*sind(t) ; ...
     -cosd(t)*sind(t) -sind(t)^2       cosd(t)*sind(t)  sind(t)^2])     ;

% Local Stiffness Matrix
% k  1    2     3    4      5    6     7    8     9      10   11
k = [k(0) k(60) k(0) k(120) k(0) k(60) k(0) k(60) k(120) k(0) k(120)];

% Element Nodes
% elem  1        2        3        4        5         6
elem = [1 2 3 4; 5 6 3 4; 5 6 7 8; 7 8 3 4; 9 10 1 2; 11 12 9 10; ...
%       7            8          9        10         11
        11 12 13 14; 1 2 13 14; 1 2 5 6; 13 14 5 6; 13 14 9 10];

% Local to global transformation
K = zeros(14);
for i = 1:11
    K(elem(i,:), elem(i,:)) = K(elem(i,:), elem(i,:)) + k(:, (i-1)*4 + 1:i*4);
end

A = .01;
E = 200e9;
k_part = K([1 2 3 4 5 6 7 9 10 13 14], [1 2 3 4 5 6 7 9 10 13 14]);
f_part = [0; 0; 0; 0; 0; -1000e3; 0; 0; 0; 0; -1000e3];
d_part = (k_part * A * E) \ f_part;

d = [d_part(1:7); 0; d_part(8:9); 0; 0; d_part(10:11)];
f = A * E * K * d;

% calculating stress and strain
epsilon = @(t,n) ((1/L) * [-cosd(t) -sind(t) cosd(t) sind(t)] * d(elem(n,:)));

% elem    1              2               3              4
strain = [epsilon(0, 1); epsilon(60, 2); epsilon(0, 3); epsilon(120, 4); ...
%         5              6               7              8                         
          epsilon(0, 5); epsilon(60, 6); epsilon(0, 7); epsilon(60, 8); ...
%         9                10              11
          epsilon(120, 9); epsilon(0, 10); epsilon(120, 11)];

stress = E * strain;








