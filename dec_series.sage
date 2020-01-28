class dec_series(object):
    def __init__(self,N,ell):
        self.N = N
        self.rN = range(1,N+1)
        self.rL = range(1,ell+1)
        self.g = function("g")

        if ell == 1:
            xv = var("x")
            self.gx = g(xv)
            self.x=[None, xv]
            self.xvars = [xv]
        else:
            self.xvars = list(var(["x" + str(i) for i in range(1,ell+1)]))
            self.x = [None] + list(self.xvars)
            self.gx = g(*self.xvars)

        
        self.Rqx = PolynomialRing(QQ, self.xvars + ["q"])
        self.Fqx = Frac(self.Rqx)

        if N == 1:
            self.gvars = [var("g")]
            self.Rg = PolynomialRing(self.Fqx,"g")
            self.qiffy_poly = self.qiffy_poly_N1
        else:
            self.gvars = var(["g" + str(i) for i in self.rN])
            self.Rg = PolynomialRing(self.Fqx,self.gvars)
        self.g = [None] + list(self.gvars)
        
        
        self.Fg = Frac(self.Rg)

        self.e = SymmetricFunctions(self.Fqx).elementary()
        
        self._R = PolynomialRing(self.Fqx, "G")
        G = self._R.gen()
        q = self.Fqx("q")
        self.Weqx0 = (1-G)^N*q - G^N

        q = var("q")
        self.Weq = self.gx^N * prod((1-xi*self.gx for xi in self.xvars))  - q*(1-self.gx)^N
  
        
    def qiffy_diff(self,f,dlist):
        fcur = f
        for xj in dlist:
            fcur = sum(( diff(fcur, gi) * self.dgdx(gi,xj) for gi in self.gvars )) + diff(fcur,xj)
        fcur0 = fcur.subs({xi:0 for xi in self.xvars})
        return self.qiffy_frac(fcur0)
    
    @cached_method
    def dgdx_(self,xj):
        return solve(diff(self.Weq,xj) == 0, diff(self.gx,xj))[0].rhs()
    
    @cached_method
    def dgdx(self,gi,xj):
        return self.dgdx_(xj).subs({self.gx:gi})

    
    def qiffy_frac(self, f):
        return self.qiffy_poly(f.numerator())/self.qiffy_poly(f.denominator())    

    def qiffy_poly(self,f):
        fsym = self.e.from_polynomial(self.Rg(f)).restrict_parts(self.N)
        result = 0
        for part, coef in fsym.monomial_coefficients().items():
            part_result = 1
            for pi in part:
                part_result *= self.Weqx0[pi-1]/self.Weqx0[self.N]
            result += coef*part_result
            
        return factor(result)

    def qiffy_poly_N1(self,f):
        q = self.Rqx("q")
        return factor(self.Rg(f).subs({self.Rg.gen(): -1/(1+q)}))