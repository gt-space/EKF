function dq = quaternion_exp(v)
% exponential map for unit quaternions
vnorm = norm(v);
if vnorm < 1e-6
    dq = [1;v];
    dq = dq/norm(dq);
else
    dq = [cos(vnorm);
          sin(vnorm)*(v/vnorm)];
end