%%
%-----------------------------------------------------
% Input to  the Model

E1 = 135e3;                                  % Youngs Modulus in 1 Direction (MPa)
E2 = 10e3;                                   % Youngs Modulus in 2 Direction (MPa)
nu12 = 0.25;                                 % Poission's Ratio
G12 = 8e3;                                   % Shear Modulus (MPa)
thick = [0.125 0.125 0.125 0.125];           % Thickness of each Ply from Bottom to Top (mm)
ply_angles = [0 45 45 0];                    % Angle of each Ply from Bottom to Top (deg)

%-----------------------------------------------------
%%

total_thick = sum(thick);                    
n = max(size(ply_angles));                   

z = zeros(1,n+1);
for i = 1:n+1
    if i == 1
        z(i) = -total_thick/2;
    else
    z(i) = z(i-1) + thick(i-1);
    end
end

nu21 = nu12*(E2/E1);
R = [1 0 0;0 1 0;0 0 2];
Q = [E1/(1-nu12*nu21) nu12*E2/(1-nu12*nu21) 0;nu12*E2/(1-nu12*nu21) E2/(1-nu12*nu21) 0;0 0 G12];

A = zeros(3,3);B = zeros(3,3);D = zeros(3,3);
for k = 1:n
    theta = ply_angles(k);
    T = [cosd(theta)^2 sind(theta)^2 2*cosd(theta)*sind(theta);sind(theta)^2 cosd(theta)^2 ...
    -2*cosd(theta)*sind(theta); -cosd(theta)*sind(theta) cosd(theta)*sind(theta) cosd(theta)^2-sind(theta)^2];
    A_temp = (inv(T)*Q*R*T*inv(R))*(z(k+1)-z(k));
    B_temp = (1/2) * (inv(T)*Q*R*T*inv(R))*(z(k+1)^2-z(k)^2);
    D_temp = (1/3) * (inv(T)*Q*R*T*inv(R))*(z(k+1)^3-z(k)^3);
    A = A+A_temp;
    B = B + B_temp;
    D = D + D_temp;
end 

ABD = [A B;B D]; % ABD Matrix
for i = 1:6
    for j = 1:6
        if abs(ABD(i,j)) < 1e-6
            ABD(i,j) = 0.0;
        end
    end
end