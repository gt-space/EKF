function CB2I = quaternion2DCM(q)
% quaternion describes rotation from inertial frame to body frame
% so DCM should convert body coordinates to inertial coordinates
% since frame and coordinates vary contravariantly

% Ensure input is a 4-element vector
if numel(q) ~= 4
    error('Input quaternion must be a 4-element vector.');
end
% Normalize the quaternion (to prevent numerical drift)
q = q / norm(q);

s = q(1);
v = q(2:4);
CB2I = (s^2 - dot(v,v))*eye(3) + 2*(v*v') + 2*s*skew(v);
% CI2B = (s^2 - dot(v,v))*eye(3) + 2*(v*v') - 2*s*skew(v);
end




