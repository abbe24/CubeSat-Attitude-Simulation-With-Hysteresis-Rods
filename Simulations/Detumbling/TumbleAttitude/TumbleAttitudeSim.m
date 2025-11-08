%{

QUESTIONS: Euler angle rates vs. Angular Velocities?, Body Fixed Frame vs
intertial Frame? 

The purpose of this simulation is to simulate a tumbling 1.5U cubesat for 100
orbits starting on AUG 1 2026 00:00:00 based on initial angular velocities using the Newton-Euler rotational 
Equations of motion for rigid-bodies. 

INPUTS: ISS Longitude Latitude Altitude ORBIT DATA (from STK), WMM Magnetic flux density data
(From STK), Moments of inertia of cubesat (IXX, IYY, IZZ), Initial Angular
Velocities (Wx, Wy, Wz), Initial Euler Angles (Phi, Theta, Psi).

Outputs: Euler angles over 100 orbits (Phi, Theta, Psi), Euler angle rates
over 100 orbits (PhiDot, ThetaDot, PsiDot), Angular Velocities in ECEF
Frame over 100 orbits

%}