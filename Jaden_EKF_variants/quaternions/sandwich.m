function y = sandwich(q,x)
% compute the sandwich product y = qxq*
    y = qmult(q,qmult(x,qconj(q)));
end