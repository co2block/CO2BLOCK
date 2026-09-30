function FD = FD_Nor(x,R,R_ext,csi, gamma)
	if x<=R		
		if R <= 0.8216*R_ext 
			FD = log(R/x);

		else
			FD =  log(R_ext/x) + 2/2.25*(R/R_ext)^2 -3/4 + 1/(2*R_ext^2)*(x^2-csi^2/gamma);				
		end
			
	else
		FD = 0;
	end	
	
end
