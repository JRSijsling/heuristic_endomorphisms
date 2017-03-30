attach("DiscriminantBound.sage")
attach("TwistPolynomials.sage")
attach("EndomorphismRankBound.sage")
attach("NonIsogenous.sage")
attach("Genus2Factors.sage")

attach("MagmaInterface.m")

def testHigherGenus() : 
    R.<u> = PolynomialRing(QQ)

    # RM by Q(\sqrt{17})
    f = R([-3, 8, 5, 7, 2, 1]) 
    d = 2*f.disc()
    C = HyperellipticCurve(f)

    # Stupid full CM
    # f = u^5 + 1
    # f = u^7 + 1
    # d = 2*f.disc()
    # C = HyperellipticCurve(f)

    # Reducible example
    # f = u^6-7*u^4+14*u^2-7
    # d = 2*f.disc()
    # C = HyperellipticCurve(f)


    # Trivial endomorphisms in genus 3
    # C = HyperellipticCurve(u^7-2*u^6+5*u^5-2*u^4+u^3+4*u^2+2*u,u^2+1)
    # d = 9449
    
    # RM by Q(sqrt(5))
    # C = HyperellipticCurve(R([0, -2, 1]), R([1, 1, 0, 1]))
    # d = 191

    LPolys = [ 0 for i in range(0,maxP) ]

    for p in range(2,maxP) :
        if is_prime(p) and d%p != 0 :
            Cp = C.base_extend(FiniteField(p))
            LPolys[p] = Cp.frobenius_polynomial()
            # print LPolys[p]
            
            
    type, bd = DiscriminantBound(LPolys, d)


    print "Type", type
    print "Bound", bd.factor()

    print "Bound on the Z-rank of the endomorphism ring", EndomorphismRankBound(LPolys, d, C.genus())
    

def TestEC() : 
    LPolys1 = [ 0 for i in range(0,maxP) ]
    LPolys2 = [ 0 for i in range(0,maxP) ]
    # p1 = [0, -1, 1, 0, 0]
    # p2 = [0, -1, 1, -10, -20]
    
    p1 = [0, -1, 0, 1, 0]
    p2 = [0, 1, 0, 1, 0]
    
    E1 = EllipticCurve(p1)
    E2 = EllipticCurve(p2)
    d = E1.discriminant() * E2.discriminant()
    for p in range(2,maxP) :
        if is_prime(p) and d%p != 0 :
            E1p = E1.base_extend(FiniteField(p))
            E2p = E2.base_extend(FiniteField(p))
            LPolys1[p] = E1p.frobenius_polynomial()
            LPolys2[p] = E2p.frobenius_polynomial()
    print "Finished computing L-functions"
    print "Can prove there is no isogeny over the ground field?", CertifyNonIsogenous(LPolys1, LPolys2, false)
    print "Can prove there is no isogeny geometrically?", CertifyNonIsogenous(LPolys1, LPolys2)
    
    
def TestPicard() :
    R.<x> = PolynomialRing(QQ)
    LPolys = magma.computeLPolys2(100)
    print "Finished computing L-polynomials"
    LPolysInternal = [R(l) for l in LPolys]
    LPolysInternal.insert(0,0)
    print "Finished converting L-polys to Sage"
    type, bound = DiscriminantBound(LPolysInternal, 1)
    
    print "Type", type
    print "Discriminant bound", bound
    
def TestDetectGenus2Factor() :
    R.<x> = PolynomialRing(QQ)

    f = x^5-x^4+x^3
    h = x^4+1
    d = 7744
    # OK, confirms that the quotient abelian surface has no extra endomorphisms

    f = x^7+x^6+x^5+x^3+x^2+x
    h = x^4+x^2+1
    d = 3993
    # Cannot prove anything, which is correct because the quotient abelian surface is QM

    
    C = HyperellipticCurve(f,h)
    LPolys = [ 0 for i in range(0,maxP) ]

    
    for p in range(2,maxP) :
        if is_prime(p) and d%p != 0 :
            Cp = C.base_extend(FiniteField(p))
            LPolys[p] = Cp.frobenius_polynomial()
            # print LPolys[p]
    
    Genus2LPolys = Genus2FactorTwistedLPolys(LPolys)
    type, bound = DiscriminantBound(Genus2LPolys, 1, true)
    print "Type", type
    print "Discriminant bound", bound
