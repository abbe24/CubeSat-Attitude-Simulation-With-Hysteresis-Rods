%{ 

The purpose of this simulation is to calculate the earths magnetic flux
density H based on eq. 17, 18, and 19 in the PMAC paper by Gerhardt. this
simulation assumes a circular orbit. 

%}

clear; clc; close all;

%-------------------------------------------------------------------------
%define constants
%-------------------------------------------------------------------------
    %permiability of free space (Mu(sub0))
FREESPACEPERMIABILITY = ((4)*pi) * (10^(-7));

%-------------------------------------------------------------------------
%defined orbital parameters
%-------------------------------------------------------------------------

    %Earth's average radius in km
Re = 6371;
    %Altitude (based on ISS avg altitude in km)
alt = 400;
    %semi-major axis of orbit in km
a = Re + alt;

    %inclination in degrees (based on ISS inclination)
i_deg = 55.0; 
    %initial argument of latitude in degrees
u0_deg = 0.0;

    %Convert inclination and initial argument of latitude to radians
i_rad = deg2rad(i_deg);
u0_rad = deg2rad(u0_deg);

%-------------------------------------------------------------------------
%Calculated magnetic parameters
%-------------------------------------------------------------------------

    %{
    Gravitational effect of earth (mu = GM) (Km^3/s^2) G is newtons constant of
    gravitation and m is mass of earth. (source - pg 198, Space Mission
    Engineering, Wertz, Everett, and puschell)
    %}
MuEarth = 398600.4418;

%{
Calculate Heq Magnetic field strength at earths equator
    Earths horizontal field intensity in Teslas (T) based on NOAA data
    (8/1/2026),(LAT 0, LONG 0, ALT 400km)(HorizontalIntensity bc mag
    field is horizontal at equator
%}
    Beq = 0.0000226068;
        %calculate Magnetic Flux density based on relationship to magnetic
        %field strength [ Heq = Beq / Mu(sub0) ]
    Heq = (Beq) / (FREESPACEPERMIABILITY);

%-------------------------------------------------------------------------
%Calculated Orbital Parameters
%-------------------------------------------------------------------------

    %calculate mean motion (movement of satellite through u in rad/s) based on (NASA Fundamentals of Orbital Mechanics
    %Ch 7 eq (7-47)) 
n = sqrt((MuEarth)/(a^3));
  
    %calculate orbital period (s)
T = (2 * pi) / n;

%-------------------------------------------------------------------------
%Simulate one orbit
%-------------------------------------------------------------------------  

    %Number of time steps
N = 1000;
    %linspace(linearly spaced vector), linspace(x1,x2,N) creats a vector of
    %N evenly spaced points between x1 and x2. Takes 1000 time steps
    %between time 0 and end of an orbit
t = linspace(0, T, N);
    %simulate argument of latitude for N(1000) evenly spaced time steps in
    %orbit. initial argument of latitude + (mean motion * time vector (1000 elements))
u = u0_rad + n .* t;

%Eqs (17),(18), and (19) from Gerhardt
H1 = 3 * Heq * sin(i_rad) * cos(i_rad) .* (sin(u).^2);
H2 = -3 * Heq * sin(i_rad) .* sin(u) .* cos(u);
H3 = Heq * (1 - 3 * (sin(i_rad).^2) .* (sin(u).^2));
