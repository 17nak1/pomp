current_params= c(R0=3.132490e+01 , amplitude=3.883620e-01 , gamma=7.305000e+01 , mu=6.469830e-04 , sigma=4.566000e+01 ,rho= 4.598709e-01 ,psi= 1.462546e-01 ,S_0= 3.399189e-02 ,E_0=2.336327e-04 ,R_0=9.657741e-01,I_0=4.221789e-07)


pfilter(m1,params=current_params,Np=100,filter.mean = T,pred.mean=T,pred.var=T,save.states = TRUE, max.fail=3000) -> ss;ss@loglik