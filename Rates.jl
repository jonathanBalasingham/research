
R_AB(alpha, beta, gamma, K, T) = alpha*(T/(300*K))^beta*exp(-gamma/T)
R_CRProton(alpha, zeta) = alpha*zeta
R_CRPhoton(alpha, beta, zeta, omega, T, K, E) = alpha*(T/(300*K)^beta *((E*zeta)/(1-omega)))
R_UV(alpha, F_UV, k, A_v) = alpha*F_UV*exp(-k*A_v)

struct Parameters
    zeta::Float64
    omega::Float64
    K::Float64
    T::Float64
    F_UV::Float64
    A_v::Float64
    k::Float64
    E::Float64
end

function calculateRates!(rdata, parameters)
    rdata[!, "rate"] .= 0.0
    for row in eachrow(rdata)
        if row.re2 == "CRP"
            row.rate = R_CRProton(row.alpha,parameters.zeta)
        elseif row.re2 == "CRPHOT"
            # E parameter is missing!!
            row.rate = R_CRPhoton(row.alpha, row.beta, parameters.zeta, parameters.omega, parameters.T,parameters.K,0.5)
        elseif row.re2 == "PHOTON"
            row.rate = R_UV(row.alpha, parameters.F_UV, parameters.k, parameters.A_v)
        else
            row.rate = R_AB(row.alpha,row.beta,row.gamma,parameters.K,parameters.T)
        end
    end
end


