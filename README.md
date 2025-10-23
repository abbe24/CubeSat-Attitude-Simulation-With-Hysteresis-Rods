# CubeSat-Attitude-Simulation-with-Hysteresis-Rods
This is a Matlab simulation project in progress for the Utah State GASRATS satellite team. The goal of GASRATS is to test a prototype of a patch antenna, which is a special antenna overlaid on a solar panel to save space. The team has opted to use hysteresis rods for passive detumbling, and myself and the other members of the ADCS subteam have produced this preliminary simulation. 

The Matlab simulation file is named "cubesat_detumble5". This simulation projects a 1U cubesat given the initial conditions to lose at least 90% of its rotation velocity by 500 sec. The simulation was made using Euler’s rigid body rotational kinematics and dynamics equations. A gimball lock on the kinematic equations is required and implemented to prevent division by zero. A photo of the used kinematic and dynamic equations is provided for reference.

Going forward, the team plans to refine the calculation for the CubeSat's damping constant and to integrate orbital mechanics and position into the simulation to, for example, provide more accurate magnetic field strength data for the position of the satellite around Earth's orbit instead of using a generalized number (like 40 microteslas).

The simulation was written with the help of Google's Gemini AI.
