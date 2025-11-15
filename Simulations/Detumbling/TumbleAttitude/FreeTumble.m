function dydt = FreeTumble(t, y, I)
% Computes the derivative of angular velocity for free tumbling (no torques)
% y = [wx; wy; wz]

wx = y(1);
wy = y(2);
wz = y(3);

Ixx = I(1);
Iyy = I(2);
Izz = I(3);

% Euler rotational equations (no torques)
wx_dot = ((Iyy - Izz)/Ixx) * wy * wz;
wy_dot = ((Izz - Ixx)/Iyy) * wz * wx;
wz_dot = ((Ixx - Iyy)/Izz) * wx * wy;

dydt = [wx_dot; wy_dot; wz_dot];
end