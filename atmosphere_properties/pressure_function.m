% uses a reference pressure p_0 and a 15th order best-fit polynomial (stored as an array of 16 coefficients) to estimate atmosphere properties at some altitude
% 0th order coefficient first, 15th order last
% turns an altitude, h (meters), into a pressure, p
function p = pressure_function(h, p_0, a)
    p = p_0 * 10^(a * [1.0; h; h^2; h^3; h^4; h^5; h^6; h^7; h^8; h^9; h^10; h^11; h^12; h^13; h^14; h^15]);
    p = p * 10^-5; % convert pascal to bar
end


% example p_0 and a values with 0.999926 r^2 over -2 km to 160 km:
% simulated for 35.3478813171387, -117.806831359863 (lat, long) at 2025-04-25 14:00:00 UTC 
% p_0 = 102715.47296217596;
% a = [-0.004366921037228173;
%    -5.0964993627948606 * 10^-5;
%    1.6319216411827967 * 10^-9;
%    -6.053135304413484 * 10^-13;
%    6.911049045104908 * 10^-17;
%    -4.516216131249696 * 10^-21;
%    1.90623861182507 * 10^-25;
%    -5.457995055935324 * 10^-30;
%    1.0906415917824724 * 10^-34;
%    -1.5456197431929747 * 10^-39;
%    1.5609417634668778 * 10^-44;
%    -1.1149608634620734 * 10^-49;
%    5.500087730582754 * 10^-55;
%    -1.7817554862391996 * 10^-60;
%    3.409963106133886 * 10^-66;
%    -2.9213421378235885 * 10^-72]


