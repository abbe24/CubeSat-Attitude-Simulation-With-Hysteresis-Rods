
clear; clc; close all;

%The purpose of this code is to explore the possibility of integrating STK
%with matlab.

%Launch STK application
app = actxserver('STK12.application');
root = app.Personality2; %DONT KNOW WHAT THIS DOES
%Make STK Application visible
app.Visible = 1;

%create Scenario from 08/01/2026 - 08/02/2026 named Matlab_Integration
scenario = root.Children.New('eScenario','MATLAB_Integration');
scenario.SetTimePeriod('1 Aug 2026 00:00:00.000','2 Aug 2026 00:00:00.000');
scenario.StartTime = '1 Aug 2026 00:00:00.000';
scenario.StopTime = '2 Aug 2026 00:00:00.000';
root.ExecuteCommand('Animate * Reset'); %DONT KNOW WHAT THIS DOES

%create satellite
satellite = scenario.Children.New('eSatellite','GASRATS');

