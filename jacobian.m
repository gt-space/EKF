clear
syms q [4 1] real
syms magI [3 1] real

h = @(q) q2v(qmult(qmult([q(1);-q(2:4)],[0;magI]),q));

Df = jacobian([h(q)],[q]) 
% Df =

[2*magI1*q1 + 2*magI2*q4 - 2*magI3*q3, 2*magI1*q2 + 2*magI2*q3 + 2*magI3*q4, 2*magI2*q2 - 2*magI1*q3 - 2*magI3*q1, 2*magI2*q1 - 2*magI1*q4 + 2*magI3*q2]
[2*magI2*q1 - 2*magI1*q4 + 2*magI3*q2, 2*magI1*q3 - 2*magI2*q2 + 2*magI3*q1, 2*magI1*q2 + 2*magI2*q3 + 2*magI3*q4, 2*magI3*q3 - 2*magI2*q4 - 2*magI1*q1]
[2*magI1*q3 - 2*magI2*q2 + 2*magI3*q1, 2*magI1*q4 - 2*magI2*q1 - 2*magI3*q2, 2*magI1*q1 + 2*magI2*q4 - 2*magI3*q3, 2*magI1*q2 + 2*magI2*q3 + 2*magI3*q4]
 