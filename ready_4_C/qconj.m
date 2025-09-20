function qbar = qconj(q)
% only for unit q
if abs(norm(q) - 1) > 1e-6
   warning('Quaternion is not unit norm. Conjugate is not inverse.')
end
qbar = [q(1); -q(2:4)];
end