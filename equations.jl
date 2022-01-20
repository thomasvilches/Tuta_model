
    #### RUnge_kutta
    function eqODE(du,u,P,t)
        du[1] = P.n*P.phi*(1-u[1]/P.K)*u[3]-P.gammaT*u[1]-P.rho*u[5]*u[1]-P.alphaE*u[6]*u[1]-P.muET*u[1] ##ET
        du[2] = P.gammaT*u[1]-P.alphaN*u[6]*u[2]-P.eta*u[2]-P.muNT*u[2] ##NT
        du[3] = P.eta*u[2]-P.muAT*u[3] ##AT
        du[4] = P.rho*u[5]*u[1]-P.gammaP*u[4]-P.alphaP*u[6]*u[4]-P.muEP*u[4] ##EP
        du[5] = P.gammaP*u[4]-P.muAP*u[5] ##AP
        du[6] = P.betaE*u[6]*u[1]+P.betaN*u[6]*u[2]+P.betaP*u[6]*u[4]-P.muAM*u[6] ##AM
    end