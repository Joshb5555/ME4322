% six bar linkage
% static equilibrium
clc
clear
% define the joints
A = [7 4 0];
B = [5 16 0];
C = [25 25 0];
D = [23 10 0];
E = [18 35 0];
F = [43 32 0];
G = [45 17 0];
% Define the lengths of the bars
lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D - C);
lBE = norm(E - B);
lEF = norm(F - E);
lFG = norm(G - F);
% weight
WAB = [0 -1 0];
WBEC = [0 -1 0];
WCD = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];
% COM of each link
% BEC is a ternary link carrying three joints, so its centroid is the
% average of all three points and divides by 3. The binary links average
% their two endpoints.
S1 = (A+B)/2;
S2 = (B+C+E)/3;
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;
syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin
ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
InputTorque = [0 0 Tin];
% Applied force
% This load acts on link FG. It is taken at S5, so it enters the force
% balance for FG but produces no moment about that link's COM. Move it and
% eqn10 needs a cross(rApplied-S5, AppliedForce) term.
AppliedForce = [50 0 0];
% Sign convention
% ForceX is the force the downstream link applies to the upstream link.
% By Newton's third law the downstream link feels -ForceX. Whichever sign
% a link uses in its force balance has to be reused in its moment balance,
% otherwise the two equations describe different free bodies and the
% assembled system comes out rank deficient with no solution.
% Static equibrium conditions for Link AB
% Sum of forces = 0
% Fa + Fb + WEightofAB = 0
eqn1 = ForceA + ForceB + WAB == 0;
% sum of moments = 0 with respect to COM of link AB
% S1A x FA +S1B x FB + InputTorque = 0
eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + InputTorque == 0;
% equations for link BEC
% -ForceB is the reaction coming back from AB at B. ForceC and ForceE are
% what CD and EF push onto this link at C and E.
eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;
% sum of moments = 0 with respect to CONM of link BEC
eqn4 = cross(B-S2,-ForceB) + cross(C-S2,ForceC) + cross(E-S2,ForceE) == 0;
% equations for link CD
eqn5 = -ForceC + ForceD + WCD == 0;
% sum of moments = 0 with repect to COM of link CD
eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;
% equations for link EF
eqn7 = -ForceE + ForceF + WEF == 0;
% sum of moments = 0 with respect to COM of link EF
eqn8 = cross(E-S4,-ForceE) + cross(F-S4,ForceF) == 0;
% equations for link FG
eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;
% sum of moments = 0 with respect to COM of link FG
% ForceF keeps the same minus sign it carries in eqn9
eqn10 = cross(F-S5,-ForceF) + cross(G-S5,ForceG) == 0;
% solving
% Each force equation contributes two useful scalars in x and y. Each
% moment equation contributes one in z. That is 15 equations against the
% 14 joint force components plus Tin, so the system is square.
eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];
StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin ]);
Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);
ForceB_x = double(StaticSolution.FBx);
ForceB_y = double(StaticSolution.FBy);
ForceC_x = double(StaticSolution.FCx);
ForceC_y = double(StaticSolution.FCy);
ForceD_x = double(StaticSolution.FDx);
ForceD_y = double(StaticSolution.FDy);
ForceE_x = double(StaticSolution.FEx);
ForceE_y = double(StaticSolution.FEy);
ForceF_x = double(StaticSolution.FFx);
ForceF_y = double(StaticSolution.FFy);
ForceG_x = double(StaticSolution.FGx);
ForceG_y = double(StaticSolution.FGy);
InputTorque_value = double(StaticSolution.Tin);
% Display the results
disp('Static Equilibrium Forces:');
disp(['Force A_x: ', num2str(Force_Ax)]);
disp(['Force A_y: ', num2str(Force_Ay)]);
disp(['Force B_x: ', num2str(ForceB_x)]);
disp(['Force B_y: ', num2str(ForceB_y)]);
disp(['Force C_x: ', num2str(ForceC_x)]);
disp(['Force C_y: ', num2str(ForceC_y)]);
disp(['Force D_x: ', num2str(ForceD_x)]);
disp(['Force D_y: ', num2str(ForceD_y)]);
disp(['Force E_x: ', num2str(ForceE_x)]);
disp(['Force E_y: ', num2str(ForceE_y)]);
disp(['Force F_x: ', num2str(ForceF_x)]);
disp(['Force F_y: ', num2str(ForceF_y)]);
disp(['Force G_x: ', num2str(ForceG_x)]);
disp(['Force G_y: ', num2str(ForceG_y)]);
disp(['Input Torque: ', num2str(InputTorque_value)]);


syms wBEC wCD
omega_AB = [0 0 1];
omega_BEC = [0 0 wBEC];
omega_CD = [0 0 wCD];

eqn11 = cross(omega_AB,B-A) + cross(omega_BEC,C-B) + cross(omega_CD,D-C) == 0;

loop1Solution = solve(eqn11, [wBEC wCD]);

% Extract angular velocities from the solution
wBEC_value = double(loop1Solution.wBEC)
wCD_value = double(loop1Solution.wCD)

omegaBEC = [0 0 wBEC_value];
omegaCD = [0 0 wCD_value];

% second loop
%DCEFGD

syms wEF wFG
omega_EF = [0 0 wEF];
omega_FG = [0 0 wFG];

% Define the angular velocities for the second loop
eqn12 = cross(omegaCD,C-D) + cross(omegaBEC,E-C) + cross(omega_EF,F-E) + cross(omega_FG,G-F) == 0;

% Solve for angular velocities in the second loop
loop2Solution = solve(eqn12, [wEF wFG]);

% Extract angular velocities from the second loop solution
wEF_value = double(loop2Solution.wEF)
wFG_value = double(loop2Solution.wFG)

% Angular acceleration
% Loop 1 ABCDA

syms aBEC aCD
% Angular acceleration for the first loop
alpha_AB = [0 0 0];
alpha_BEC = [0 0 aBEC];
alpha_CD = [0 0 aCD];

a_B_A = cross(alpha_AB,B-A) + cross(omega_AB,cross(omega_AB,B-A));
a_C_B = cross(alpha_BEC,C-B) + cross(omegaBEC,cross(omegaBEC,C-B));
a_D_C = cross(alpha_CD,D-C) + cross(omegaCD,cross(omegaCD,D-C));

eqn13 = a_B_A + a_C_B + a_D_C == 0;

loop1AccSolution = solve(eqn13,[aBEC aCD]);

alphaBEC = double(loop1AccSolution.aBEC)
alphaCD = double(loop1AccSolution.aCD)

alphaBEC_vector = [0 0 alphaBEC];
alphaCD_vector = [0 0 alphaCD];

syms aEF aFG

alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

a_C_D = cross(alphaCD_vector,C-D) + cross(omegaCD,cross(omegaCD,C-D));
a_E_C = cross(alphaBEC_vector,E-C) + cross(omegaBEC,cross(omegaBEC,E-C));

ang_Vel_EF = [0 0 angularVelocity_EF];
ang_vel_FG = [0 0 angularVelocity_FG];

a_F_E = cross(alpha_EF,F-E) + cross()
