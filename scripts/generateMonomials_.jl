# as_termes oder as_coifichens_dict




using SymPy

@syms X1 X2     # define symbolic variables
sqrt=√X1        # overwrite square root function with symbolic square root


# polynomial order
grad = [1,3,6,10,15]
order = 2

# basis functions phiᵢ
phi = [
  sqrt(2);
  2 - 6*X1;
  2*sqrt(3)*(1 - X1 - 2*X2);
  sqrt(6)*((10*X1 - 8).*X1 + 1);
  sqrt(3)*((5*X1 - 4).*X1 + (-15*X2 + 12).*X2 - 1);
  3*sqrt(5)*((3*X1 + 8*X2 - 4).*X1 + (3*X2 - 4).*X2 + 1);
  2*sqrt(2)*(-1+(15+(-45+35*X1).*X1).*X1);
  2*sqrt(6)*(-1+(13+(-33+21*X1).*X1).*X1+(2+(-24+42*X1).*X1).*X2);
  2*sqrt(10)*(-1+(9+(-15+7*X1).*X1).*X1+(6+(-48+42*X1).*X1+(-6+42*X1).*X2).*X2);
  2*sqrt(14)*(-1+(3+(-3+X1).*X1).*X1+(12+(-24+12*X1).*X1+(-30+30*X1+20*X2).*X2).*X2);
  sqrt(10)*(1+(-24+(126+(-224+126*X1).*X1).*X1).*X1);
  sqrt(30)*(1+(-22+(105+(-168+84*X1).*X1).*X1).*X1+(-2+(42+(-168+168*X1).*X1).*X1).*X2);
  5*sqrt(2)*(1+(-18+(69+(-88+36*X1).*X1).*X1).*X1+(-6+(102+(-312+216*X1).*X1).*X1+(6+(-96+216*X1).*X1).*X2).*X2);
  sqrt(70)*(1+(-12+(30+(-28+9*X1).*X1).*X1).*X1+(-12+(132+(-228+108*X1).*X1).*X1+(30+(-300+270*X1).*X1+(-20+180*X1).*X2).*X2).*X2);
  3*sqrt(10)*(1+(-4+(6+(-4+X1).*X1).*X1).*X1+(-20+(60+(-60+20*X1).*X1).*X1+(90+(-180+90*X1).*X1+(-140+140*X1+70*X2).*X2).*X2).*X2);
]

d₁phi = diff(phi,X1)    # first partial derivative of phiᵢ
d₂phi = diff(phi,X2)    # second partial derivative of phiᵢ

monomials=Set{String}()

for i in 1:grad[order]
    for j in 1:grad[order]
        for k in 1:grad[order]
            pol = expand(phi[i]*phi[j]*d₁phi[k])
# print("Polynom: ")
#println(pol)

            polString = split(replace(replace(string(pol),"- ","-"),"+ ",""))
            regex = r"^(?(?=.*?\*X)(?'coefficient'.*?)\*(?'monomial'X?1?\^?\d?\*?X?2?\^?\d?)|(?'const'.*))$"

            monsMatch = map(x->match(regex,x), polString)
            for m in monsMatch
              if nothing == m["const"]
                  mon = m["monomial"]
                  coeff = m["coefficient"]
              else
                  mon = "1"
                  coeff = m["const"]
              end
                push!(monomials, mon)
# println("monomials: " * mon)
            end

            #TODO
            # pol2 = Poly(phi[i]*phi[j]*d₂phi[k])

        end
    end
end

println("The following monomials appear when using polynomials up to order $(order) for 'dphi*phi*phi':")
for i in monomials
    println(i)
end
