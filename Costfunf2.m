function f = Costfunf2(newsol,volfrac,S1,S2,S)
    [intrusion_max,A1]=LS_dyn_FEm(newsol,S1,S2,S);
     penalfact = 200000;
    f = intrusion_max+penalfact*max((A1/volfrac)-1,0);
end

