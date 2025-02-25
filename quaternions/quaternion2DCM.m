function CB2I = quaternion2DCM(q)
% quaternion describes rotation from inertial frame to body frame
% so DCM should convert body coordinates to inertial coordinates
% since frame and coordinates vary contravariantly
s = q(1);
v = q(2:4);
CB2I = (s^2 - dot(v,v))*eye(3) + 2*(v*v') + 2*s*skew(v);
% CI2B = (s^2 - dot(v,v))*eye(3) + 2*(v*v') - 2*s*skew(v);
end




