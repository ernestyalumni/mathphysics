import sympy as sp
# 4D hyperspherical coords on S^3: k = |k| (cos a, sin a cos b, sin a sin b cos c, sin a sin b sin c)
a,b,c = sp.symbols('a b c', real=True)
n0 = sp.cos(a); n1 = sp.sin(a)*sp.cos(b)
dOmega = sp.sin(a)**2*sp.sin(b)
def avg(f):
    num = sp.integrate(sp.integrate(sp.integrate(f*dOmega,(c,0,2*sp.pi)),(b,0,sp.pi)),(a,0,sp.pi))
    return sp.simplify(num/(2*sp.pi**2))
area = sp.integrate(sp.integrate(sp.integrate(dOmega,(c,0,2*sp.pi)),(b,0,sp.pi)),(a,0,sp.pi))
print("area S^3 =", area)
A1 = avg(n1**2); A2 = avg(n1**2*n0**2)
print("<n1^2> =", A1, " <n1^2 n0^2> =", A2)
C = sp.symbols('C', positive=True)
h_avg = 4*C*(A1 - 4*A2 - 1)
print("angular average of h |k|^4 =", sp.simplify(h_avg))
L = sp.symbols('L', positive=True); k = sp.symbols('k', positive=True)
shell = sp.integrate(area*k**3/k**4,(k,1/L,1))
print("shell integral of |k|^-4 =", sp.simplify(shell))
kappa = sp.simplify(-h_avg*shell/(2*sp.pi)**4/sp.log(L))
print("kappa_G =", kappa, " expected 11*C/(24*pi^2):", sp.simplify(kappa-11*C/(24*sp.pi**2))==0)
# Casimir constant for su(N) in the proof's normalization <X,Y> = -(1/N) Re tr(XY), tr(adX adY) = -C <X,Y>
for N in (2,3):
    # basis of su(N): anti-Hermitian
    import itertools
    basis=[]
    for i in range(N):
        for j in range(i+1,N):
            E=sp.zeros(N); E[i,j]=1; E[j,i]=-1; basis.append(E)
            F=sp.zeros(N); F[i,j]=sp.I; F[j,i]=sp.I; basis.append(F)
    for i in range(N-1):
        D=sp.zeros(N); D[i,i]=sp.I; D[i+1,i+1]=-sp.I; basis.append(D)
    dim=len(basis)
    G=sp.Matrix(dim,dim,lambda i,j: -sp.Rational(1,N)*sp.re((basis[i]*basis[j]).trace()))
    def coords(M):
        v=sp.Matrix([-sp.Rational(1,N)*sp.re((basis[i]*M).trace()) for i in range(dim)])
        return G.LUsolve(v)
    ad=[sp.Matrix.hstack(*[coords(X*Y-Y*X) for Y in basis]) for X in basis]
    K=sp.Matrix(dim,dim,lambda i,j:(ad[i]*ad[j]).trace())
    ratio=sp.simplify(-K[0,0]/G[0,0])
    print(f"su({N}): dim={dim}, C_g={ratio}, Killing=-C*G everywhere:", sp.simplify(K+ratio*G)==sp.zeros(dim), " expected 2N^2 =",2*N*N)
