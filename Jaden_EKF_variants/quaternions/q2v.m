function v = q2v(q)
if abs(q(1)) > 1e-6
   warning('q2v should only be used on pure quaternions.')
end
v = q(2:4);
end