function qbar = qconj(q)
% only for unit q
% if abs(norm(q) - 1) > 1e-6
%    error('Error. quaternion is not unit norm.\n')
% end
qbar = [q(1); -q(2:4)];
end