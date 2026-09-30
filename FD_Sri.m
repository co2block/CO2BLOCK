function FD = FD_Sri(r,R,R_ext,csi, gamma)
u = @(x,Xc) 2.25/4*(x/Xc)^2;
Sr_app = @(G) (log(0.5614594835668851698/G)+0.9713*G - 0.1742*G^2); 
%Srivastava approx for well function replaces 2*log(R/r); here the algorithm has been modified to be adapted to closed reservoir
	if u(r,R)<=1 
		if R <= 0.8216*R_ext
			FD = 1/2*Sr_app(u(r,R));
		else         
			FD =  log(R_ext/r) + 2/2.25*(R/R_ext)^2 -3/4 + 1/(2*R_ext^2)*(r^2-csi^2/gamma);
		end
	else
		FD = 1/2*(1/(u(r,R)*exp(u(r,R)))* (u(r,R)+0.3637)/(u(r,R)+1.282));
	end
end

% Reference
%Srivastava, R., & Guzman, A. (1998). Practical Approximations of the Well Function. Groundwater, 36(5).
% u = (r/(4/3*R))^2;
% u_top =  (r_top/(4/3*R))^2;
% u_bot =  (r_bot/(4/3*R))^2;
% Sr_app = (log(0.5614594835668851698/u)+0.9713*u - 0.1742*u^2); 
% Sr_app_top = (log(0.5614594835668851698/u_top)+0.9713*u_top - 0.1742*u_top^2); 
% Sr_app_bot = (log(0.5614594835668851698/u_bot)+0.9713*u_bot - 0.1742*u_bot^2); 
% Srivastava approx for well function replaces 2*log(R/r)