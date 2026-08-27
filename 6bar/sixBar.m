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

% Pack the loop 2 results the same way loop 1 was packed. Everything below
% uses these two names. The earlier draft reached for angularVelocity_EF
% and angularVelocity_FG, which were never assigned, so it errored here.
omegaEF = [0 0 wEF_value];
omegaFG = [0 0 wFG_value];

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

% Angular acceleration
% Loop 2 DCEFGD, the same circuit that produced wEF and wFG

syms aEF aFG

alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

% The first two terms are already numeric because loop 1 has been solved.
% Only a_F_E and a_G_F still carry unknowns.
a_C_D = cross(alphaCD_vector,C-D) + cross(omegaCD,cross(omegaCD,C-D));
a_E_C = cross(alphaBEC_vector,E-C) + cross(omegaBEC,cross(omegaBEC,E-C));
a_F_E = cross(alpha_EF,F-E) + cross(omegaEF,cross(omegaEF,F-E));
a_G_F = cross(alpha_FG,G-F) + cross(omegaFG,cross(omegaFG,G-F));

eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == 0;

loop2AccSolution = solve(eqn14,[aEF aFG]);

alphaEF = double(loop2AccSolution.aEF)
alphaFG = double(loop2AccSolution.aFG)

alphaEF_vector = [0 0 alphaEF];
alphaFG_vector = [0 0 alphaFG];

alphaAB_vector = [0 0 0];

disp(' ');
disp('Angular Velocities (rad/s, positive counterclockwise):');
disp(['omega AB:  ', num2str(1)]);
disp(['omega BEC: ', num2str(wBEC_value)]);
disp(['omega CD:  ', num2str(wCD_value)]);
disp(['omega EF:  ', num2str(wEF_value)]);
disp(['omega FG:  ', num2str(wFG_value)]);

disp(' ');
disp('Angular Accelerations (rad/s^2):');
disp(['alpha AB:  ', num2str(0)]);
disp(['alpha BEC: ', num2str(alphaBEC)]);
disp(['alpha CD:  ', num2str(alphaCD)]);
disp(['alpha EF:  ', num2str(alphaEF)]);
disp(['alpha FG:  ', num2str(alphaFG)]);

% Velocities at the joints
% A, D and G are pinned to ground so they do not move. Every other point
% is reached with v = omega x r, where r runs from a point of known
% velocity on the same link to the point of interest. Any point can be
% reached by more than one route because the loops close, and the two
% routes have to agree. That agreement is the check on the loop solutions.
velocityA = [0 0 0];
velocityD = [0 0 0];
velocityG = [0 0 0];

% B belongs to AB, which pivots about A
velocityB = cross(omega_AB, B-A);
% C belongs to CD, which pivots about D
velocityC = cross(omegaCD, C-D);
% E belongs to BEC. B is the point on that link whose velocity is known
velocityE = velocityB + cross(omegaBEC, E-B);
% F belongs to FG, which pivots about G
velocityF = cross(omegaFG, F-G);

% Cross checks along the other branch. Both residuals should be zero to
% machine precision. If they are not, loop 1 or loop 2 was mis-assembled.
checkC = velocityB + cross(omegaBEC, C-B) - velocityC;
checkF = velocityE + cross(omegaEF, F-E) - velocityF;

disp(' ');
disp('Joint Velocities [vx vy vz]:');
disp(['Velocity A: ', num2str(velocityA)]);
disp(['Velocity B: ', num2str(velocityB)]);
disp(['Velocity C: ', num2str(velocityC)]);
disp(['Velocity D: ', num2str(velocityD)]);
disp(['Velocity E: ', num2str(velocityE)]);
disp(['Velocity F: ', num2str(velocityF)]);
disp(['Velocity G: ', num2str(velocityG)]);
disp(['Closure residual at C: ', num2str(norm(checkC))]);
disp(['Closure residual at F: ', num2str(norm(checkF))]);

% Accelerations at the joints
% a = alpha x r + omega x (omega x r). The second term is the centripetal
% part and always points back toward the reference point.
accelerationA = [0 0 0];
accelerationD = [0 0 0];
accelerationG = [0 0 0];

accelerationB = cross(alphaAB_vector, B-A) + cross(omega_AB, cross(omega_AB, B-A));
accelerationC = cross(alphaCD_vector, C-D) + cross(omegaCD, cross(omegaCD, C-D));
accelerationE = accelerationB + cross(alphaBEC_vector, E-B) + cross(omegaBEC, cross(omegaBEC, E-B));
accelerationF = cross(alphaFG_vector, F-G) + cross(omegaFG, cross(omegaFG, F-G));

checkAccC = accelerationB + cross(alphaBEC_vector, C-B) + cross(omegaBEC, cross(omegaBEC, C-B)) - accelerationC;
checkAccF = accelerationE + cross(alphaEF_vector, F-E) + cross(omegaEF, cross(omegaEF, F-E)) - accelerationF;

disp(' ');
disp('Joint Accelerations [ax ay az]:');
disp(['Acceleration A: ', num2str(accelerationA)]);
disp(['Acceleration B: ', num2str(accelerationB)]);
disp(['Acceleration C: ', num2str(accelerationC)]);
disp(['Acceleration D: ', num2str(accelerationD)]);
disp(['Acceleration E: ', num2str(accelerationE)]);
disp(['Acceleration F: ', num2str(accelerationF)]);
disp(['Acceleration G: ', num2str(accelerationG)]);
disp(['Closure residual at C: ', num2str(norm(checkAccC))]);
disp(['Closure residual at F: ', num2str(norm(checkAccF))]);

% Velocities at the centers of mass
% Same rule as the joints. For a link that pivots on ground the reference
% is the ground pin, for a floating link it is whichever joint already has
% a velocity. S2 sits on BEC and is referenced to B, S4 sits on EF and is
% referenced to E.
velocityS1 = cross(omega_AB, S1-A);
velocityS2 = velocityB + cross(omegaBEC, S2-B);
velocityS3 = cross(omegaCD, S3-D);
velocityS4 = velocityE + cross(omegaEF, S4-E);
velocityS5 = cross(omegaFG, S5-G);

disp(' ');
disp('Center of Mass Velocities [vx vy vz]:');
disp(['Velocity S1 (AB):  ', num2str(velocityS1)]);
disp(['Velocity S2 (BEC): ', num2str(velocityS2)]);
disp(['Velocity S3 (CD):  ', num2str(velocityS3)]);
disp(['Velocity S4 (EF):  ', num2str(velocityS4)]);
disp(['Velocity S5 (FG):  ', num2str(velocityS5)]);

% Accelerations at the centers of mass
accelerationS1 = cross(alphaAB_vector, S1-A) + cross(omega_AB, cross(omega_AB, S1-A));
accelerationS2 = accelerationB + cross(alphaBEC_vector, S2-B) + cross(omegaBEC, cross(omegaBEC, S2-B));
accelerationS3 = cross(alphaCD_vector, S3-D) + cross(omegaCD, cross(omegaCD, S3-D));
accelerationS4 = accelerationE + cross(alphaEF_vector, S4-E) + cross(omegaEF, cross(omegaEF, S4-E));
accelerationS5 = cross(alphaFG_vector, S5-G) + cross(omegaFG, cross(omegaFG, S5-G));

disp(' ');
disp('Center of Mass Accelerations [ax ay az]:');
disp(['Acceleration S1 (AB):  ', num2str(accelerationS1)]);
disp(['Acceleration S2 (BEC): ', num2str(accelerationS2)]);
disp(['Acceleration S3 (CD):  ', num2str(accelerationS3)]);
disp(['Acceleration S4 (EF):  ', num2str(accelerationS4)]);
disp(['Acceleration S5 (FG):  ', num2str(accelerationS5)]);