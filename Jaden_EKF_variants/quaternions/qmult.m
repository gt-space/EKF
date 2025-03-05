function q_out = qmult(q1, q2)
    % Quaternion multiplication: q_out = q1 * q2
    % Assumes quaternions are in [scalar, vector] form: [q0, q1, q2, q3]
    
    % Extract components
    s1 = q1(1); v1 = q1(2:4);  % Scalar and vector part of q1
    s2 = q2(1); v2 = q2(2:4);  % Scalar and vector part of q2
    
    % Compute quaternion product
    s_out = s1 * s2 - dot(v1, v2);  % Scalar part
    v_out = s1 * v2 + s2 * v1 + cross(v1, v2);  % Vector part
    
    % Combine into output quaternion
    q_out = [s_out; v_out];
end