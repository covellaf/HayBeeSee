%% Initialisation script for jumping robot Simscape model
clear
clf

Gravity = 9.81;
AirDensity = 1.225;

% Mass Properties
BodyMass = 1.0;
LandingLegMass = 0.1;
JumpingLegMass = 0.3;

% Contact Forces
LandingLegDynamicFriction = 0.5;
LandingLegStaticFriction  = 0.5;
JumpingLegDynamicFriction = 2;
JumpingLegStaticFriction  = 2;
FloorPlaneWidth  = 1;
FloorPlaneLength = 10;
FloorPlaneOffset = 2;

% Dimensions
ArmLength = 1.5;
ArmWidth  = 0.1;
LandingLegHeight = 0.5;
JumpingLegLength = 1.2;
JumpingLegCompressedLength = JumpingLegLength-0.35;
JumpingLegAttachPoint = -0.05;%1;
JumpingLegStartAngle = rad2deg(asin(LandingLegHeight/JumpingLegCompressedLength));
% JumpingLegStartAngle = rad2deg(asin((LandingLegHeight-0.05)/JumpingLegCompressedLength));
% JumpingBendAngle     = rad2deg(asin(0.5*JumpingLegCompressedLength/(0.5*JumpingLegLength)));
% JumpingLegStartAngle = -(90-JumpingBendAngle-JumpingLegStartAngle);
HardStop = 170;
% JumpingLegStartAngle = 5;

HipHingeStiffness = 0.02;
HipHingeDamping = 0.01;
Bar.Length0             = 0.65;
Bar.MaxCompressionRatio = (Bar.Length0-(JumpingLegLength-JumpingLegCompressedLength))/Bar.Length0;
Bar.MinCompressionRatio = 1;


Bar.Radius      = 1.5e-3;
Bar.Width       = 24e-3;
Bar.E           = 140e9;%112.4e9;
% Bar.I           = ((Bar.Radius)^3)*Bar.Width/12;
Bar.I           = Bar.Radius^4*pi/4;
Bar.Density     = 1600;
Bar.Yield       = 2.05e9;
Bar.Number      = 2;
Bar.Bundle      = 6;
Bar.Euler       = Bar.Number*Bar.Bundle*pi^2*Bar.E*Bar.I/Bar.Length0^2;
P_Eu            = Bar.Number*Bar.Bundle*pi^2*Bar.E*Bar.I/Bar.Length0^2;
Bar.Resolution  = 20;
Bar.CompressionSteps    = 10;

theta0 = sin(linspace(0,3*pi/2,Bar.Resolution))';
A = zeros(3,Bar.Resolution);
A(1,1) = 0;
A(2,Bar.Resolution) = 0;
b = [0; 0; 0];

Bar.Theta       = zeros(Bar.Resolution,Bar.CompressionSteps);
Bar.Y           = zeros(Bar.Resolution,Bar.CompressionSteps);
Bar.X           = zeros(Bar.Resolution,Bar.CompressionSteps);
Bar.Alpha0      = zeros(1,Bar.CompressionSteps);
Bar.YMax        = zeros(1,Bar.CompressionSteps);
Bar.MaxStress   = zeros(1,Bar.CompressionSteps);
Bar.Energy      = zeros(1,Bar.CompressionSteps);
Bar.Length      = linspace(Bar.Length0,...
                           Bar.MaxCompressionRatio*Bar.Length0,...
                           Bar.CompressionSteps);

for j = Bar.CompressionSteps:-1:1

[Bar.Theta(:,j),Bar.Energy(j)] = ...
                fmincon(@(thet) elasten(thet,Bar.Length0,Bar.Radius,...
                Bar.E,Bar.Resolution,Bar.I),theta0,[],[],A,b,[],[],@(thet)...
                geomcon(thet,Bar.Length0,Bar.Resolution,...
                0,Bar.Length(j),0));
            
theta0 = Bar.Theta(:,j);

Bar.X(:,j) = cumsum((Bar.Length0/Bar.Resolution)*cos(Bar.Theta(:,j)));
Bar.Y(:,j) = cumsum((Bar.Length0/Bar.Resolution)*sin(Bar.Theta(:,j)));
Bar.YMax(:,j)   = max(Bar.Y(:,j));                
Bar.Energy(j)   = Bar.Energy(j)*Bar.Number*Bar.Bundle;
Bar.Alpha0(j)   = Bar.Theta(Bar.Resolution,j);
Bar.MaxStress   = max(Bar.Radius*Bar.E*abs(diff(Bar.Theta))...
                    /(Bar.Length0/Bar.Resolution));

                
end

Bar.Force = abs(mean(diff(Bar.Energy(2:end))./diff(Bar.Length(2:end))));
max(Bar.MaxStress/Bar.Yield)
max(Bar.Energy)

subplot(1,2,1)
plot(Bar.X(:,Bar.CompressionSteps),Bar.Y(:,Bar.CompressionSteps),'ko-')
hold on
xlabel('X (m)')
ylabel('Y (m)')
axis equal
grid on

load('JumpLegTest.mat')


subplot(1,2,2)
plot(Bar.Length,[0 abs(diff(Bar.Energy)./diff(Bar.Length))])
hold on
plot([Bar.Length(1) Bar.Length(end)],[P_Eu P_Eu])
plot([Bar.Length(1) Bar.Length(end)],[Bar.Force Bar.Force])

plot(JumpLegTest(:,2)*1e-3+max(Bar.Length)+0.014,-JumpLegTest(:,3)*1e3)

xlabel('Length (m)')
ylabel('Force (N)')
grid on
shg

P_Eu  =Bar.Force;

%% Other Design Definitions

% Elastica Leg Definition
Beam.E   = 120e9;
Beam.rho = 1.6e3;
Beam.G   = 10e9;
Beam.W   = 4*12e-3;
Beam.T   = 2*2e-3;
Beam.I   = (Beam.W*Beam.T^3)/12;
Beam.L   = 2;

% BucklingLeg
% w = 3e-3;
% b = 50e-3;
% I = w^3*b/12;
% E = 120e9;
% L = JumpingLegLength;
% nbars = 2;
% P_Eu = nbars*pi^2*E*I/L^2;


% Motor
MotorVoltage       = 16;
MotorStallTorque   = 2.32;
MotorFreeSpeed     = 2027;
ArmatureResistance = 3.89;
GearRatio          = 30;
SEA_Spring         = 1e6;
PEASpring          = 1.2*360/pi;
PEASpringEquilib   = 30;


% Four Bar Leg Initialisation
Scale = 2;
Knee = 5*Scale;
Shin = 30*Scale;
Base = 5*Scale;
BaseAngle   = deg2rad(30);
BaseY = Base*sin(BaseAngle);
BaseX = Base*cos(BaseAngle);
FrontThigh  = 30*Scale;
RearThigh   = 29.5*Scale;
InitialAngle = -113;
