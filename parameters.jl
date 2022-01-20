


@with_kw mutable struct EDOparameters @deftype Float64

    n = 0.55
    phi = 14.0
    K = 400
    gammaT = 1/4.0
    muET = 0.045
    muAT = 1/16.0
    muNT = 0.0316
    eta = 1/16.0
    rho = 0.0985
    gammaP = 1/10.0
    muEP = 0.00913
    muAP = 1/13.0
    p_t = 0.01
    alphaE = 19.47
    alphaN = 0.15##not sure
    alphaP = 0.001 #not sure
    muAM = 1/35
    p_pa = 0.00125


    betaE = p_t*alphaE
    betaN = p_t*alphaN
    betaP = p_pa*alphaP

    step::Int64 = 100

end