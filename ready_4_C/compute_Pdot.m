function Pdot = compute_Pdot(x,P,Q,a_meas,w_meas,we)

F = compute_F(x,a_meas,w_meas,we);
G = compute_G(x);

Pdot = F*P + P*F' + G*Q*G';

end