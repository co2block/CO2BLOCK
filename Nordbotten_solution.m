
function [PD] = Nordbotten_solution(r,R,psi,R_ext,gamma)
    
if r <= psi

    PD = gamma*log(psi/r) + FD_Nor(psi,R,R_ext);
      
elseif (r > psi) && (r <= R)
    PD = FD_Nor(r,R,R_ext);

else

    PD = 0;
end
end