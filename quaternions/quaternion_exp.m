function dq = quaternion_exp(v)
% exponential map for unit quaternions
    vnorm = norm(v);
    dq = [cos(vnorm);
          (sin(vnorm)/vnorm)*v];
end