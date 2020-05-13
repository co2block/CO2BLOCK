function FD = FD_Nor(x,R,R_ext)
if x < R_ext
    
    if R <= R_ext
        FD = log(R/x);

    else
        FD =  log(R_ext/x) + 2/2.25*(R/R_ext)^2 -3/4 ;
    end
    
else
    FD = 0;
        
end
end