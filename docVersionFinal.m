%% Setup

clc; clear;

% Scale parameters
h = 1E-1;

x = [0, 2];
% x = [0.95, 1.05];
% x = [1.000, 1.002];
xPts = linspace(x(1), x(end), 500);

y = [-3*h, 0];
% y = [-10*h, 0];

% initialize theoretical comparisons
yTheory = 0;
yTheory2 = 0;
yTheory3 = 0;

% choose N by running exactly one of the three boxes below

%% N = 2

% Parameters
l = 6; C1 = 1; C2 = 1; b1 = 2.0; b2 = 2.0;

% Algebraic equation for z
zEqn = @(z) (2*i.*z - C1*h^b1).*(2*i.*z - C2*h^b2).*exp(-2*i*l.*z./h) ...
    - C1*C2*h^(b1+b2);

% Comparison to theory
yTheory = ((b1+b2)/(2*l))*h*log(h) + h/(2*l)*log(abs(C1*C2) ./ (4.*xPts.^2));

%% N = 3

% Parameters
C1 = 1; C2 = 1; C3 = 1;
% l1 = 4; l2 = 2; b1 = 0.5; b2 = 2.0; b3 = 2.0;
l1 = 4; l2 = 2; b1 = 0.5; b2 = 0.5; b3 = 2.0;

% Algebraic equation for z
zEqn = @(z) (2*i.*z - C1*h^b1).*(2*i.*z - C2*h^b2).*(2*i.*z .... 
    - C3*h^b3).*exp(-2*i*(l1+l2).*z./h) ...
    - (C1*h^b1).*(C2*h^b2).*(2*i.*z - C3*h^b3).*exp(-2*i*(l2).*z./h) ...
    - (2*i.*z - C1*h^b1).*(C2*h^b2).*(C3*h^b3).*exp(-2*i*(l1).*z./h) ...
    - (C1*h^b1).*(2*i.*z + C2*h^b2).*(C3*h^b3);

% Comparison to theory
% yTheory = (((b1+b3)/2)/(l1+l2))*h*log(h) ...
%     + h/(2*(l1+l2))*log(abs(C1*C3) ./ (4.*xPts.^2)); % gamma 13
yTheory = ((b1+b2)/(2*l1))*h*log(h) ...
    + h/(2*l1)*log(abs(C1*C2) ./ (4.*xPts.^2)); % gamma 12
yTheory2 = ((b3-b2)/(2*l2))*h*log(h) ...
    + h/(2*l2)*log(abs(C3/C2)) .* ones(size(xPts)); % gamma 23

%% N = 4

% Parameters
C1 = 1; C2 = 1; C3 = 1; C4 = 1;
% l1 = 2; l2 = 3; l3 = 1; b1 = 2; b2 = 2; b3 = 2; b4 = 2; % gamma 14
% l1 = 2; l2 = 3; l3 = 1; b1 = 2; b2 = 2; b3 = .5; b4 = 2; % gamma 13, 34
l1 = 2; l2 = 3; l3 = 1; b1 = 2; b2 = .5; b3 = .5; b4 = 2; % gamma 23, 12, 34
% l1 = 3; l2 = 1; l3 = 2; b1 = .5; b2 = .5; b3 = 2; b4 = 3; % gamma 12, 24
% l1 = 3; l2 = 2; l3 = 1; b1 = .5; b2 = .5; b3 = 2; b4 = 3; % gamma 12, 23, 34
% l1 = 1; l2 = 3; l3 = 1; b1 = 2; b2 = .5; b3 = .5; b4 = 2; % double string
% l1 = 2; l2 = 3; l3 = 1; b1 = 1.5; b2 = .5; b3 = .5; b4 = 2; % crossed strings
% C1 = 10; C3 = -5; l1 = 2; l2 = 3; l3 = 1; b1 = 2; b2 = .5; b3 = .5; b4 = 2;

% Algebraic equation for z
% For convenience, the scaling is different than for N=2,3
R1 = @(z) (C1*h^b1)./(2*i.*z-C1*h^b1);
R2 = @(z) (C1*h^b2)./(2*i.*z-C2*h^b2);
R3 = @(z) (C1*h^b3)./(2*i.*z-C3*h^b3);
R4 = @(z) (C1*h^b4)./(2*i.*z-C4*h^b4);
zEqn = @(z) -exp(-2*i*(l1+l2+l3).*z./h) ...
    + R1(z).*R2(z).*exp(-2*i*(l2+l3).*z./h) ...
    + R2(z).*R3(z).*exp(-2*i*(l1+l3).*z./h) ... 
    + R3(z).*R4(z).*exp(-2*i*(l1+l2).*z./h) ...
    + R2(z).*(1+2*R3(z)).*R4(z).*exp(-2*i*(l1).*z./h) ... 
    - R1(z).*R2(z).*R3(z).*R4(z).*exp(-2*i*(l2).*z./h) ...
    + R1(z).*(1+2*R2(z)).*R3(z).*exp(-2*i*(l3).*z./h) ...
    + R1(z).*(1+2*R2(z)).*(1+2*R3(z)).*R4(z);

% Comparison to theory

% Dominant strings
% yTheory = ((b1+b2)/(2*l1))*h*log(h) ...
%     + h/(2*l1)*log(abs(C1*C2) ./ (4.*xPts.^2)); % gamma 12
yTheory = ((b2+b3)/(2*l2))*h*log(h) ...
    + h/(2*l2)*log(abs(C2*C3) ./ (4.*xPts.^2)); % gamma 23
% yTheory = ((b1+b3)/(2*(l1+l2)))*h*log(h) ...
%     + h/(2*(l1+l2))*log(abs(C1*C3)./ (4.*xPts.^2)); % gamma 13
% yTheory = ((b1+b4)/(2*(l1+l2+l3)))*h*log(h) ...
%     + h/(2*(l1+l2+l3))*log(abs(C1*C4) ./ (4.*xPts.^2)); % gamma 14

% Nondominant strings
% yTheory2 = ((b4-b2)/(2*(l2+l3)))*h*log(h) * ones(size(xPts)); % gamma 24
yTheory2 = ((b4-b3)/(2*l3))*h*log(h) * ones(size(xPts)); % gamma 34
% yTheory3 = ((b3-b2)/(2*(l2)))*h*log(h) * ones(size(xPts)); % gamma 23
yTheory3 = ((b1-b2)/(2*l1))*h*log(h) * ones(size(xPts)); % gamma 21

%% Compute contours and intersections

zReal = @(x, y) real(zEqn(x + i*y));
zImag = @(x, y) imag(zEqn(x + i*y));

% Compute 0 level contours of the real and imaginary part
mesh = 2500;
figure()
f = fcontour(zReal, [x(1) x(end) y(1) y(end)], 'LevelList', 0, 'MeshDensity', mesh);
Mreal = f.ContourLines.VertexData;
close;
figure()
f = fcontour(zImag, [x(1) x(end) y(1) y(end)], 'LevelList', 0, 'MeshDensity', mesh);
Mimag = f.ContourLines.VertexData;
close;

% Find the intersections of the above contours
% uses intersections() from Fast and Robust Curve Intersections by Douglas Schwarz
Mreal = Mreal(:,1:3:end); % uniformly remove data for speed
Mimag = Mimag(:,1:3:end);
epsilon = 0.03*h; % cutoff intersections within epsilon of domain
Mimag(1, abs(Mimag(1,:) - x(1)) < epsilon) = NaN; 
Mimag(1, abs(Mimag(1,:) - x(end)) < epsilon) = NaN;
Mimag(2, abs(Mimag(2,:) - y(1)) < epsilon) = NaN; 
Mimag(2, abs(Mimag(2,:) - y(end)) < epsilon) = NaN;
Mreal(1, abs(Mreal(1,:) - x(1)) < epsilon) = NaN; 
Mreal(1, abs(Mreal(1,:) - x(end)) < epsilon) = NaN;
Mreal(2, abs(Mreal(2,:) - y(1)) < epsilon) = NaN; 
Mreal(2, abs(Mreal(2,:) - y(end)) < epsilon) = NaN;
[resReal, resImag, iout, jout] = intersections(Mreal(1,:), Mreal(2,:), ...
    Mimag(1,:), Mimag(2,:), 0);

% Diagnostics: check the green dots are at the intersection of contours
figure()
hold on;
plot(Mreal(1, :), Mreal(2, :), '.', 'Color', 'red')
plot(Mimag(1, :), Mimag(2, :), '.', 'Color', 'blue')
plot(resReal, resImag, '.', 'Color', 'green', 'MarkerSize', 20)
xlim(x);
ylim(y);

% Diagnostics: check a resonance is the root with << h^2 precision of the main equation
resIndex = 15;
fprintf("f(%.16f + %.16f*i) = %.1E + %.1E i\n", resReal(resIndex), resImag(resIndex), ...
    zReal(resReal(resIndex), resImag(resIndex)), zImag(resReal(resIndex), resImag(resIndex)))

%% Plotting

figure()
hold on; box on; pbaspect([1 1 1]);
set(gca, 'DefaultLineLineWidth', 2, 'DefaultTextInterpreter', 'latex', ...
    'fontsize', 16);
plot(resReal, resImag, 'x', 'Color', 'red', 'MarkerSize', 10)
plot(xPts, yTheory, 'Color', 'blue');
plot(xPts, yTheory2, 'Color', 'blue');
plot(xPts, yTheory3, 'Color', 'blue');

xlim(x);
ylim([-0.3,0]);
xlabel("Re z");
ylabel("Im z");
legend("resonances $z$", "theory", "location", "southeast", ...
    "Interpreter", "latex", 'fontsize', 20)
