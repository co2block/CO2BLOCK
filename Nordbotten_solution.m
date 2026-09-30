function [PD] = Nordbotten_solution(r,R,csi,L,gamma, Sriv)
delta = gamma/(1-gamma);
r_bot = csi*sqrt(gamma);
r_top = csi/sqrt(gamma);

if Sriv == 'on'
	FD = @(a,b,c,d,e) FD_Sri(a,b,c,d,e);     
else
	FD = @(a,b,c,d,e) FD_Nor(a,b,c,d,e);  
end

if r > R_ext
	PD = 0
else
	if R>csi
		if r <= r_bot
			PD = gamma*log(r_bot/r)   + sqrt(gamma)/csi*(r_top-r_bot) +  FD(r_top,R,R_ext,csi,gamma);
			
		elseif (r > r_bot) && (r <= r_top)
			PD = sqrt(gamma)/csi*(r_top-r) + FD(r_top,R,R_ext,csi,gamma); 
			
		elseif (r > r_top) 
			PD = FD(r,R,R_ext,csi,gamma);
				
		end
	else
		PD = gamma* FD(r,R,R_ext,csi,gamma); 
end
end
end


