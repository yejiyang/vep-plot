function prin_stress_example

sigmax=12000; sigmay=15000; tauxy=8000;
original_stress=[sigmax tauxy; tauxy sigmay];

[Vectors, Principal] = eig (original_stress); %Eigen value function in MATLAB

principal_stress = Principal %Principal stresses
Value = Vectors (2,1)/Vectors (1,1); % Slope of eigenvector is a ratio of its (y/x)-components
RadianOrientation = atan(Value);
PrincipalOrientation = RadianOrientation * (180/pi) %Principal stress orientation in degrees

ShearTheta = RadianOrientation - (pi/4);
MaxShearOrientation = ShearTheta*(180/pi) %Maximum shear stress orientation in degrees

T=[cos(ShearTheta) sin(ShearTheta); sin(ShearTheta) cos(ShearTheta)]; % Transformation matrix

Maxshear = T*original_stress*transpose(T) % Calculating Maximum shear using Matrix method
