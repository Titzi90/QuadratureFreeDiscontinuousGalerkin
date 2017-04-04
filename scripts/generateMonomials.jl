# as_termes oder as_coifichens_dict




using SymPy

@syms X1 X2     # define symbolic variables


# polynomial order
grad = [1,3,6,10,15]
order = 2

# basis functions phiᵢ
phi = [
    sqrt(Sym(2));                                                                                                                                        # phi_1:  c          -> P_0
    2 - 6*X1;                                                                                                                                            # phi_2:  X1
    2*sqrt(Sym(3))*(1 - X1 - 2*X2);                                                                                                                      # phi_3:  X2         -> P_1
    sqrt(Sym(6))*((10*X1 - 8).*X1 + 1);                                                                                                                  # phi_4:  X1^2
    sqrt(Sym(3))*((5*X1 - 4).*X1 + (-15*X2 + 12).*X2 - 1);                                                                                               # phi_5:  X2^2
    3*sqrt(Sym(5))*((3*X1 + 8*X2 - 4).*X1 + (3*X2 - 4).*X2 + 1);                                                                                         # phi_6:  X1*X2      -> P_2
    2*sqrt(Sym(2))*(-1+(15+(-45+35*X1).*X1).*X1);                                                                                                        # phi_7:  X1^3
    2*sqrt(Sym(6))*(-1+(13+(-33+21*X1).*X1).*X1+(2+(-24+42*X1).*X1).*X2);                                                                                # phi_8:  X1^2*X2
    2*sqrt(Sym(10))*(-1+(9+(-15+7*X1).*X1).*X1+(6+(-48+42*X1).*X1+(-6+42*X1).*X2).*X2);                                                                  # phi_9:  X1*X"^2"
    2*sqrt(Sym(14))*(-1+(3+(-3+X1).*X1).*X1+(12+(-24+12*X1).*X1+(-30+30*X1+20*X2).*X2).*X2);                                                             # phi_10: X2^3       -> P_3
    sqrt(Sym(10))*(1+(-24+(126+(-224+126*X1).*X1).*X1).*X1);                                                                                             # phi_11: x1^4
    sqrt(Sym(30))*(1+(-22+(105+(-168+84*X1).*X1).*X1).*X1+(-2+(42+(-168+168*X1).*X1).*X1).*X2);                                                          # phi_12: X1^3*X2
    5*sqrt(Sym(2))*(1+(-18+(69+(-88+36*X1).*X1).*X1).*X1+(-6+(102+(-312+216*X1).*X1).*X1+(6+(-96+216*X1).*X1).*X2).*X2);                                 # phi_13: X1^2*X2^2
    sqrt(Sym(70))*(1+(-12+(30+(-28+9*X1).*X1).*X1).*X1+(-12+(132+(-228+108*X1).*X1).*X1+(30+(-300+270*X1).*X1+(-20+180*X1).*X2).*X2).*X2);               # phi_14: X1*X2^3
    3*sqrt(Sym(10))*(1+(-4+(6+(-4+X1).*X1).*X1).*X1+(-20+(60+(-60+20*X1).*X1).*X1+(90+(-180+90*X1).*X1+(-140+140*X1+70*X2).*X2).*X2).*X2);               # phi_15: X2^4       -> P_4
]




""" calculates ...
    phis are basis functions of the polynomial space
    order is the order of polynomial space spanned by phis
"""
function intDPhiPhiPhi(phis, order)
    dPhis = (diff(phis,X1), diff(phis,X2) )   # partial derivatives of phi

    monomials=Set{SymPy.Sym}()
    orderPol = max(0, order-1+order+order)
    M = (orderPol+1)*(orderPol+2)/2
    N = length(phis)    # order of polynomial space
    @assert N == (order+1)*(order+2)/2
    #G = Array{Array{SymPy.Sym,1}}(2, N, N, N)
    #G = Array{SymPy.Sym}}(2, N, N, N, M)
    G = Array{Dict{SymPy.Sym,SymPy.Sym}}(2, N, N, N)
    monomials=Set{SymPy.Sym}()


    for m in 1:2
        dₘPhis = dPhis[m]
        for i in 1:N
            dₘPhiᵢ = dₘPhis[i]
            for j in 1:N
                phiⱼ = phis[j]
                for k in 1:N
                    phiₖ = phis[k]

                    pol = expand(dₘPhiᵢ * phiⱼ * phiₖ)

                    #G[m,i,j,k]=Array{SympPy.Sym}(0)
                    G[m,i,j,k]=Dict{SymPy.Sym,SymPy.Sym}()

                    for l in as_terms(pol)[1]
                        term = l[1]
                        ex = l[2][2]
                        if 0 == length(ex)
                            ex1 = 0
                            ex2 = 0
                        elseif 1 == length(ex)
                            ex1 = ex[1]
                            ex2 = 0
                        else
                            ex1 = ex[1]
                            ex2 = ex[2]
                        end

                        monom = X1^ex1*X2^ex2
                        coeff = term/monom
                        #intOverRefTriangle = integrate(integrate(monom, X2, 0, 1-X1), X1, 0, 1)

                        #push!(G[m,i,j,k], coeff)
                        #G[m,i,j,k,    ] = coeff
                        G[m,i,j,k][monom] = coeff
                        push!(monomials, monom)
                    end

                end
            end
        end
    end


    return G, monomials
end



function computeCoefficentMatrix(coefficents::Array{Dict{SymPy.Sym,SymPy.Sym},4}, monomials::Array{SymPy.Sym,1})
    M = size(coefficents)[1]
    N = size(coefficents)[2]
    L = length(monomials)

    G = Array{SymPy.Sym}(2, N, N, N, L)

    for m in 1:M
        for i in 1:N
            for j in 1:N
                for k in 1:N
                    for l in 1:L
                        G[m,i,j,k,l] = get(coefficents[m,i,j,k], monomials[l], Sym(0))
                        # if haskey(coefficents[m,i,j,k], monomials[l])
                        #     G[m,i,j,k,l]= coefficents[m,i,j,k][monomials[l]]
                        #end
                    end
                end
            end
        end
    end

    return G, monomials
end

function computeCoefficentMatrix(foo::Tuple{Array{Dict{SymPy.Sym,SymPy.Sym},4},Set{SymPy.Sym}})
    computeCoefficentMatrix(foo[1], collect(foo[2]))
end

# integrate array of polymials ofer ref triangle
function integrateRefTriangle(monomials::Array{SymPy.Sym,1})
    return map(i -> integrate(integrate(i, X2,0,1-X1), X1,0,1),monomials)
end
