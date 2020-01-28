class dec_series(object):
    def __init__(N,ell):
        self.N = N
        self.rN = range(1,N+1)
        self.rL = range(1,ell+1)
        if ell == 1:
            xv = var("x")
            self.gx = g(xv)
            self.x=[None, xv]
            self.xvars = [xv]
        else:
            self.xvars = var(["x" + str(i) for i in range(1,ell+1)])
            self.x = None + list(xvars)
            self.gx = g(*xvars)

        self.Rqx = PolynomialRing(QQ, self.xvars + ["q"])
        self.Fqx = Frac(self.Rqx)

        self.gvars = var(["g" + str(i) for i in rN])
        self.g = [None] + list(gvars)
        
        self.Rg = PolynomialRing(Fqx,self.gvars)
        self.Fg = Frac(self.Rg)

        self.e = SymmetricFunctions(self.Fqx).elementary()
        self._R = PolynomialRing(self.Fqx, "G")
        G = self._R.gen()

        q = self.Fqx("q")
        self.Weqx0 = (1-G)^N*q - G^N

        q = var("q")
        self.Weq = gx^N * prod((1-xi*gx for xi in self.xvars))  - q*(1-gx)^N
  
        
    def qiffy_diff(self,f,dlist):
        fcur = f
        for xj in dlist:
            fcur = sum(( diff(fcur, gi) * dgdx(gi,xj) for gi in self.gvars )) + diff(fcur,xj)
        fcur0 = fcur.subs({xi:0 for xi in self.xvars})

        return qiffy_frac(Fg(fcur0))
    
    @cache_method
    def dgdx_(xj):
        return solve(diff(self.Weq,xj) == 0, diff(gx,xj))[0].rhs()
    
    @cache_method
    def dgdx(gi,xj)
        return dgdx_(xj).subs({gx:gi})

    