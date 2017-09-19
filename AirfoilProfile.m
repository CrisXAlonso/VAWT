%% Airfoil Paths and Profiles for SolidWorks

% By: Timothy A. Coda
% Date of Last Revision: 10 Sept 2017
%Modified on September 19th

% This script is used to generate text files of XYZ coordinates to be
% imported into SolidWorks for the purpose of creating CAD models of
% airfoils/blades for the VAWT project, based out of the FDRL under the
% direction of Professor Charles Williamson.

% There are multiple ways to create CAD models of a VAWT blade.  The method
% I ultimately settled upon utilizes the "Lofted Boss/Base" tool in
% SolidWorks.  At minimum this tool requires two curves, which indicate the
% limits of the extrusion.  What makes this tool useful for our purposes is
% that it allows the selection of a "Guide Curve," which is what allows for
% move complex blade design to be modeled.  Thus, this script produces text
% files that fall into one of two categories:

% (1) Airfoil chord line, used as the "Guide Curve" for a loft extrusion in
% SolidWorks.  This is done by first specifying the desired chord line in
% polar coordines (defining r, theta, and z) and then converting them to
% Cartesian coordinates that can be read in SolidWorks.  The conversion is
% done using MATLAB's pol2cart function.
% (2) The airfoil profile (i.e. NACA 0015, NACA 0021, etc.).  This is done
% using the NACA airfoil profile equation, found at:
% https://en.wikipedia.org/wiki/NACA_airfoil

% Three different blade designs are considered below:
% (1) ZigZag
% (2) Witch of Agnesi (WoA)
% (3) Sigmoid

% *** ALL UNITS IN cm ***
close all
clear all
clc
%% ZigZag Airfoil (ZZ)

% First, specify the value of r of the chord line.  This is simply the
% radius of the turbine, the radial distance from the center of the turbine 
% at which the airfoil is mounted.  For the 1x turbine, R = 8.25 cm.
R = 8.25; % recall all units are cm

% Now, specify the blade span, which will be in the z-direction.
% Historically, the span of the blade is 22 cm:
SPAN = 22;

% It is necessary to divide the blade, and chord line, into discrete segments.
% SolidWorks has a difficult time extruding along a complex paths.
% Instead, discretize the blade based upon its features.  The ZigZag blade
% lends itself to being treated as two portions that are symmetric about
% the mid-chord point.  In general, try to define paths that are
% monotonically increasing or decreasing in nature.
% Furthermore, as will be seen below, this discretization is necessary so
% that the blade curvature is modeled in SolidWorks as desired.
n = 100; % number or discrete points comprising the chord path
Z1_ZZ = linspace(0,SPAN/2,n)';
Z2_ZZ = linspace(SPAN/2,SPAN,n)';

% Now define the corresponding theta value for each segment.  For the ZigZag
% blade with a sweep angle (SA) of 45 degrees, the theta variation is simply
% linear:
SA_ZZ = pi/4; % in radians
theta1_ZZ = linspace(0,SA_ZZ,n)';
theta2_ZZ = linspace(SA_ZZ,0,n)';

% Finally, define the radial values of the chord line.  We are not
% considering a bulbous or canted blade, so r is constant:
r1_ZZ = R*ones(1,n)';
r2_ZZ = R*ones(1,n)';

% Now we can use the pol2cart function.  Note that all of the arrays above
% are column arrays; this is the required format for the pol2cart function.
% The function returns three separate column arrays of the x, y, and z
% values.  These arrays are joined into a single matrix and then exported
% as a text file, which can be read by SolidWorks

% Additionally, it is also best to shift the curve such that it begins at
% the origin (this will make it more compatible with the airfoil curves 
% generated below).  All that needs to be done is subtraction of the value 
% of R from the x-values of the curves:
[x1_ZZ,y1_ZZ,z1_ZZ] = pol2cart(theta1_ZZ,r1_ZZ,Z1_ZZ);
x1_ZZ = x1_ZZ - R;
XYZ1_ZZ = [x1_ZZ,y1_ZZ,z1_ZZ];
[x2_ZZ,y2_ZZ,z2_ZZ] = pol2cart(theta2_ZZ,r2_ZZ,Z2_ZZ);
x2_ZZ = x2_ZZ - R;
XYZ2_ZZ = [x2_ZZ,y2_ZZ,z2_ZZ];

%% Witch of Agnesi (WoA)

% Note that this is not exactly the chord line guide curves used in the WoA
% blades made in Spring 2017; those blades consisted for 4 segments, 2 of
% which were straight.  Below only two segments are considered, thus this is very
% similar to the ZigZag blade above.  The z and r components of the guide
% curves will be the same as those above; only the theta component changes.

% The WoA function:

% f(z) = (8*a^3)/(z^2 + 4*a^2)

% The nature of the WoA equation requires shifting and scaling to achieve
% the desired profile with the desired sweep angle.  The variable "a" is
% a tuning parameter of sorts, and determines the "sharpness" of the WoA
% curve.  Furthermore, the WoA curve is
% symmetric about z = 0.  However, the other curves and profiles generated
% in this code are based in the z = [0,SPAN] range, not the 
% z = [-SPAN/2,SPAN/2] range.  Thus, when using the WoA equation,
% the independent variable (z) must be shifted to get the desired theta
% variation.
Z1_WoA = Z1_ZZ; Z2_WoA = Z2_ZZ; r1_WoA = r1_ZZ; r2_WoA = r2_ZZ;
SA_WoA = pi/4;
a = pi/4; % a "tuning parameter"
% WoA function
theta1_WoA = (8*a^3)./(((Z1_WoA-(0.5*SPAN)).^2 + 4*a^2));
% First shift:
theta1_WoA = theta1_WoA - min(theta1_WoA);
% Then scale:
beta = SA_WoA/max(theta1_WoA);
theta1_WoA = beta.*theta1_WoA;
% WoA function
theta2_WoA = (8*a^3)./(((Z2_WoA-(0.5*SPAN)).^2 + 4*a^2)); 
% First shift:
theta2_WoA = theta2_WoA - min(theta2_WoA);
% Then scale:
theta2_WoA = beta.*theta2_WoA;
[x1_WoA,y1_WoA,z1_WoA] = pol2cart(theta1_WoA,r1_WoA,Z1_WoA);
x1_WoA = x1_WoA - R;
XYZ1_WoA = [x1_WoA,y1_WoA,z1_WoA];
[x2_WoA,y2_WoA,z2_WoA] = pol2cart(theta2_WoA,r2_WoA,Z2_WoA);
x2_WoA = x2_WoA - R;
XYZ2_WoA = [x2_WoA,y2_WoA,z2_WoA];

%% Sigmoid Path

% The Sigmoid blade is the most complex; it is discretized into 5 different
% segments.  However, three of those segments are simple, straight lines.
% Thus, it is not fully necessary to define and export those segments, as
% the Loft tool will not be necessary.  Instead, a simple linear extrusion
% can be done.

% The Sigmoid function:

% f(z) = 1/(1 + exp(-b*z))

% Like the WoA curve, the Sigmoid curve will require shifting and scaling to achieve the
% desired profile.  And, like the WoA function, the Sigmoid function is
% centered about z = 0.  So, the dependent variable (z) will have to be
% shifted.  "b" represents a tuning parameter and controls the "steepness"
% of the Sigmoid curve.

% Finally, the curves created below were not actually used to create the
% "Sigmoid" blade created in Spring 2017.  In fact, the Sigmoid function
% was not even used.  The theta variation was simply linear.  However, the
% code below produces a theta variation using the Sigmoid function.

% The blade has three straight segments, all with a length of 31/6 cm, and
% two Sigmoid segments, both with a length of 3.25 cm (summing of the five
% lengths gives our SPAN value of 22 cm)

SA_sig = pi/4;
b = 1.5;
delz1 = 31/6; delz2 = 3.25;
Z1_sig = linspace(0,delz1,n)';
Z2_sig = linspace(delz1,delz1+delz2,n)';
Z3_sig = linspace(delz1+delz2,delz1+delz2+delz1,n)';
Z4_sig = linspace(delz1+delz2+delz1,delz1+delz2+delz1+delz2,n)';
Z5_sig = linspace(delz1+delz2+delz1+delz2,delz1+delz2+delz1+delz2+delz1,n)';
r1_sig = R*ones(length(Z1_sig),1);
r2_sig = R*ones(length(Z2_sig),1);
r3_sig = R*ones(length(Z3_sig),1);
r4_sig = R*ones(length(Z4_sig),1);
r5_sig = R*ones(length(Z5_sig),1);

% Sigmoid Function
theta2_sig = SA_sig.*(1./(1+exp(-b*(Z2_sig - delz1 - 0.5*delz2))));
% First shift:
theta2_sig = theta2_sig - min(theta2_sig);
% Then scale:
gamma = SA_sig/max(theta2_sig);
theta2_sig = gamma*theta2_sig;
% The theta variation for the other Sigmoid segment can be easily derived
% from theta_sig2 by reversing the order of the entries of the array,
% essentially mirroring the variation:
theta4_sig = flipud(theta2_sig);

% theta is constant for the other three segments (theta = 0 for the first
% and fifth segments, and theta = pi/4 for the third segment)
theta1_sig = zeros(length(Z1_sig),1);
theta3_sig = (pi/4)*ones(length(Z3_sig),1);
theta5_sig = zeros(length(Z5_sig),1);

% Given the simplicity of the straight segments, we will not bother
% converting them into Cartesian coordinates, nor will we write them to text
% files.  We will, however, use them for plotting to check our work.

[x2_sig,y2_sig,z2_sig] = pol2cart(theta2_sig,r2_sig,Z2_sig);
x2_sig = x2_sig - R;
XYZ2_sig = [x2_sig,y2_sig,z2_sig];
[x4_sig,y4_sig,z4_sig] = pol2cart(theta4_sig,r4_sig,Z4_sig);
x4_sig = x4_sig - R;
XYZ4_sig = [x4_sig,y4_sig,z4_sig];

%% Plot Guide Curve

% For more complex paths it may be better to first check here that the path
% is defined correctly.  Thic can be done by plotting theta versus z on a
% Cartesian plane:
figure(1)
grid on
box on
hold on
axis([0 SPAN 0 SA_sig])
plot(Z1_sig,theta1_sig,'b-','LineWidth',2)
plot(Z2_sig,theta2_sig,'r-','LineWidth',2)
plot(Z3_sig,theta3_sig,'g-','LineWidth',2)
plot(Z4_sig,theta4_sig,'c-','LineWidth',2)
plot(Z5_sig,theta5_sig,'m-','LineWidth',2)
legend('Segment 1','Segment 2','Segment 3','Segment 4','Segment 5','location','best')
xlabel('Z (cm)')
ylabel('\theta (rad)')
set(gca,'FontSize',24)
hold off
%% Airfoil Profile

% Define the chord length (in cm)
c = 6; % chord length

% Define the thickness of the airfoil as a fraction (recall that the last
% two digits of the NACA 00xx series indicate the thickness of the airfoil
% as a percent of the chord length.  So, if you want a NACA 0015, define
% the thickness as 0.15).

thickness = 0.15;

% Generate "baseline/reference" airfoil profile:
t = thickness*c;
x_air = linspace(0,c,1000)';
y_air = 5*t.*(0.2969.*sqrt(x_air./c) - 0.1260.*(x_air./c) - 0.3516.*(x_air./c).^2 + 0.2843.*(x_air./c).^3 - 0.1015.*(x_air./c).^4);
x_air = [x_air;flipud(x_air(1:end))];
y_air = [y_air;-flipud(y_air(1:end))];
z_air = zeros(1,length(x_air))';
XYZ_air = [x_air,y_air,z_air];

% It is possible that you may need to rotate the coordinates of the airfoil
% profile.  I have found that rotating the baseline profile -90 degrees makes
% it most compatible with the guide curves generated.
rot = -pi/2;
R1 = [cos(rot),-sin(rot),0;sin(rot),cos(rot),0;0,0,1];
XYZ_air = R1*XYZ_air';
XYZ_air = XYZ_air';

% Now, we must generate additional airfoil profiles to be located at the end of
% our lofted paths.  This requires translation and/or rotation of the
% "reference" airfoil profile made above.  The rotational component is derived from the sweep 
% angle; the translational component is derived from the sweep angle and
% the span of the chord line path

% At most, the number of airfoil profiles needed (including the reference
% profile) should be equal to the number of segments plus 1.  However, in
% many instances you may not need this many.  Simple, linear extrustions
% may not require the use of the Loft tool.  Furthermore, leverage the
% symmetry of the blade by using the "Mirror" tool in SolidWorks. The ZZ
% and WoA blades require one additional airfoil profile, whereas the
% Sigmoid blade requires 2.

% First, let us shift the chord line path(s) such that they are aligned
% with the airfoil trailing edge (I found that this produced better models
% SolidWorks).  To achieve this, simply translate the curves in the
% negative y-direction:
XYZ1_ZZ(:,2) = XYZ1_ZZ(:,2) - c;
XYZ2_ZZ(:,2) = XYZ2_ZZ(:,2) - c;
XYZ1_WoA(:,2) = XYZ1_WoA(:,2) - c;
XYZ2_WoA(:,2) = XYZ2_WoA(:,2) - c;
XYZ2_sig(:,2) = XYZ2_sig(:,2) - c;
XYZ4_sig(:,2) = XYZ4_sig(:,2) - c;

% Now we can generate our additional airfoil profiles:
% Second airfoil profile for ZigZag and WoA blades:
XYZ_air2_ZZ = XYZ_air;
% Rotation:
R2_ZZ = [cos(SA_ZZ),-sin(SA_ZZ),0;sin(SA_ZZ),cos(SA_ZZ),0;0,0,1];
XYZ_air2_ZZ = R2_ZZ*XYZ_air2_ZZ';
XYZ_air2_ZZ = XYZ_air2_ZZ';
% Translation:
XYZ_air2_ZZ(:,1) = XYZ_air2_ZZ(:,1) - c*sin(SA_ZZ) + XYZ1_ZZ(end,1);
XYZ_air2_ZZ(:,2) = XYZ_air2_ZZ(:,2) + c*cos(SA_ZZ) + XYZ1_ZZ(end,2);
XYZ_air2_ZZ(:,3) = XYZ_air2_ZZ(:,3) + XYZ1_ZZ(end,3);

% Second airfoil profile for Sigmoid blade:
XYZ_air2_sig = XYZ_air;
% Translation:
XYZ_air2_sig(:,3) = XYZ_air2_sig(:,3) + Z5_sig(1);
% Third airfoil profile for Sigmoid blade:
XYZ_air3_sig = XYZ_air;
% Rotation:
R3_sig = [cos(SA_sig),-sin(SA_sig),0;sin(SA_sig),cos(SA_sig),0;0,0,1];
XYZ_air3_sig = R3_sig*XYZ_air3_sig';
XYZ_air3_sig = XYZ_air3_sig';
% Translation:
XYZ_air3_sig(:,1) = XYZ_air3_sig(:,1) - c*sin(SA_sig) + XYZ2_sig(end,1);
XYZ_air3_sig(:,2) = XYZ_air3_sig(:,2) + c*cos(SA_sig) + XYZ2_sig(end,2);
XYZ_air3_sig(:,3) = XYZ_air3_sig(:,3) + XYZ2_sig(end,3);

% Export airfoil coordinates to a text file for SolidWorks:
dlmwrite('XYZAirfoil.txt',XYZ_air);
dlmwrite('XYZAirfoil2_ZZ.txt',XYZ_air2_ZZ);
dlmwrite('XYZAirfoil2_sig.txt',XYZ_air2_sig);
dlmwrite('XYZAirfoil3_sig.txt',XYZ_air3_sig);

% Export chord path coordinates to text files for SolidWorks
dlmwrite('ZigZagPath1.txt',XYZ1_ZZ);
dlmwrite('ZigZagPath2.txt',XYZ2_ZZ);
dlmwrite('WoAPath1.txt',XYZ1_WoA);
dlmwrite('WoAPath2.txt',XYZ2_WoA);
dlmwrite('SigmoidPath2.txt',XYZ2_sig);
dlmwrite('SigmoidPath4.txt',XYZ4_sig);