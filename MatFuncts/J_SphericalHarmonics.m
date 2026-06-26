function [Jx,Jy,Jz] = J_SphericalHarmonics(j)
% Pauli spin operators, spin-0, spin-1, and spin-2 operators
m_all = j:-1:-j;
Jz = diag(m_all);
 
rrank = length(m_all);
 
Jplus = zeros(rrank, rrank);
Jminus = zeros(rrank, rrank);
 
for iorbital = 2:rrank
  m = m_all(iorbital);
  Jplus(iorbital-1, iorbital) = sqrt( j*(j+1) - m*(m+1) );
end
 
for iorbital = 1:rrank-1
  m = m_all(iorbital);
  Jminus(iorbital+1, iorbital) = sqrt( j*(j+1) - m*(m-1) );
end
 
Jx = (Jplus + Jminus) / 2.0;
Jy = (Jplus - Jminus) / (2.0*1i);
 
end
