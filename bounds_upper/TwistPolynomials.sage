def twistPolynomial(poly, k) :
    R.<x,v> = PolynomialRing(ZZ, 2)
    f = R(poly)
    g = v-x^k
    T.<v> = PolynomialRing(QQ)
    return T(f.resultant(g))

def twistPolynomialOld(poly, k) :
    if k == 12 :
        return twistPoly12(poly)
        
    g=Integers()(poly.degree() / 2);
    
    SymmetricRing = PolynomialRing(QQ, 'a', 2*g)
    a = SymmetricRing.gens()
    
    P=PolynomialRing(QQ, 'x', 2*g)
    x = P.gens()
    R.<t>=PolynomialRing(P)
    
    LPoly=1;
    originalCoefficients = [ poly.coefficients ( sparse = false )[i] for i in range(0,2*g) ]
    twistedCoefficients = [ 0 for i in range(0,2*g+1) ] ;

    for i in range(0,2*g) :
        LPoly = LPoly * (t-R(x[i])^k);
        
    for j in range(0,2*g) : 
        fSym = LPoly.coefficients(sparse=false)[j];
        S = SymmetricFunctions(QQ)
        SingleCoefficient = S.from_polynomial(fSym)
        Elementary = S.elementary()
        SymmetricPoly = Elementary(SingleCoefficient)


        
        TwistedjthCoefficient = sum( c*prod( a[2 * g - i] if i<=2*g else 0 for i in part ) for (part,c) in SymmetricPoly )

        print originalCoefficients
        print TwistedjthCoefficient
        
        twistedCoefficients[j] = TwistedjthCoefficient( originalCoefficients )
        
        print twistedCoefficients[j]

    twistedCoefficients[2*g]=1
    print twistedCoefficients
    return poly.parent()(twistedCoefficients);
        

        
        
def twistPoly12(charPoly) :
    s4=charPoly.coefficients(sparse=false)[0]
    s3=-charPoly.coefficients(sparse=false)[1]
    s2=charPoly.coefficients(sparse=false)[2]
    s1=-charPoly.coefficients(sparse=false)[3]
    n4=s4^12;
    
    n3 = -4* s4^9 + 3* s4^8* (s1)^4 + 36* s4^8* (s1)^2* (s2) + 18* s4^8* (s2)^2 - 24* s4^7* (s1)^2* (s2)^3 - 12* s4^7* (s2)^4 + 2* s4^6* (s2)^6 + 36* s4^8* (s1)* (s3) - 48* s4^7* (s1)^3* (s2)* (s3) - 144* s4^7* (s1)* (s2)^2* (s3) + 60* s4^6* (s1)* (s2)^4* (s3) - 72* s4^7* (s1)^2* (s3)^2 - 72* s4^7* (s2)* (s3)^2 + 180* s4^6* (s1)^2* (s2)^2* (s3)^2 + 120* s4^6* (s2)^3* (s3)^2 - 36* s4^5* (s2)^5* (s3)^2 + 40* s4^6* (s1)^3* (s3)^3 + 240* s4^6* (s1)* (s2)* (s3)^3 - 240* s4^5* (s1)* (s2)^3* (s3)^3 + 30* s4^6* (s3)^4 - 180* s4^5* (s1)^2* (s2)* (s3)^4 - 180* s4^5* (s2)^2* (s3)^4 + 105* s4^4* (s2)^4* (s3)^4 - 72* s4^5* (s1)* (s3)^5 + 252* s4^4* (s1)* (s2)^2* (s3)^5 + 42* s4^4* (s1)^2* (s3)^6 + 84* s4^4* (s2)* (s3)^6 - 112* s4^3* (s2)^3* (s3)^6 - 96* s4^3* (s1)* (s2)* (s3)^7 - 12* s4^3* (s3)^8 + 54 * s4^2* (s2)^2* (s3)^8 + 12 * s4^2* (s1)* (s3)^9 - 12* s4* (s2)* (s3)^10 + (s3)^12
   
    
    n1 = -4 * s4^3 + 30 * s4^2 * (s1)^4 - 12 * s4 * (s1)^8 + (s1)^12 - 72 * s4^2 * (s1)^2 * (s2) + 84 * s4 * (s1)^6 * (s2) - 12 * (s1)^10 * (s2) + 18 * s4^2 * (s2)^2 - 180 * s4 * (s1)^4 * (s2)^2 + 54 * (s1)^8 * (s2)^2 + 120 * s4 * (s1)^2 * (s2)^3 - 112 * (s1)^6 * (s2)^3 - 12 * s4 * (s2)^4 + 105 * (s1)^4 * (s2)^4 - 36 * (s1)^2 * (s2)^5 + 2 * (s2)^6 + 36 * s4^2 * (s1) * (s3) - 72 * s4 * (s1)^5 * (s3) + 12 * (s1)^9 * (s3) + 240 * s4 * (s1)^3 * (s2) * (s3) - 96 * (s1)^7 * (s2) * (s3) - 144 * s4 * (s1) * (s2)^2 * (s3) + 252 * (s1)^5 * (s2)^2 * (s3) - 240 * (s1)^3 * (s2)^3 * (s3) + 60 * (s1) * (s2)^4 * (s3) - 72 * s4 * (s1)^2 * (s3)^2 + 42 * (s1)^6 * (s3)^2 + 36 * s4 * (s2) * (s3)^2 - 180 * (s1)^4 * (s2) * (s3)^2 + 180 * (s1)^2 * (s2)^2 * (s3)^2 - 24 * (s2)^3 * (s3)^2 + 40 * (s1)^3 * (s3)^3 - 48 * (s1) * (s2) * (s3)^3 + 3 * (s3)^4


    n2 = 6 * s4^6 + 48 * s4^5 * (s1)^4 + 3 * s4^4 * (s1)^8 + 36 * s4^5 * (s1)^2 * (s2) + 24 * s4^4 * (s1)^6 * (s2) - 36 * s4^5 * (s2)^2 + 72 * s4^4 * (s1)^4 * (s2)^2 - 96 * s4^4 * (s1)^2 * (s2)^3 + 40 * s4^3 * (s1)^6 * (s2)^3 + 105 * s4^4 * (s2)^4 - 60 * s4^3 * (s1)^4 * (s2)^4 + 144 * s4^3 * (s1)^2 * (s2)^5 - 112 * s4^3 * (s2)^6 + 42 * s4^2 * (s1)^4 * (s2)^6 - 72 * s4^2 * (s1)^2 * (s2)^7 + 54 * s4^2 * (s2)^8 + 12 * s4 * (s1)^2 * (s2)^9 - 12 * s4 * (s2)^10 + (s2)^12 - 72 * s4^5 * (s1) * (s3) - 144 * s4^4 * (s1)^5 * (s3) - 192 * s4^4 * (s1)^3 * (s2) * (s3) - 48 * s4^3 * (s1)^7 * (s2) * (s3) - 36 * s4^4 * (s1) * (s2)^2 * (s3) - 144 * s4^3 * (s1)^5 * (s2)^2 * (s3) - 192 * s4^3 * (s1)^3 * (s2)^3 * (s3) + 96 * s4^3 * (s1) * (s2)^4 * (s3) - 180 * s4^2 * (s1)^5 * (s2)^4 * (s3) + 144 * s4^2 * (s1)^3 * (s2)^5 * (s3) - 144 * s4^2 * (s1) * (s2)^6 * (s3) - 96 * s4 * (s1)^3 * (s2)^7 * (s3) + 72 * s4 * (s1) * (s2)^8 * (s3) - 12 * (s1) * (s2)^10 * (s3) + 306 * s4^4 * (s1)^2 * (s3)^2 + 120 * s4^3 * (s1)^6 * (s3)^2 + 36 * s4^4 * (s2) * (s3)^2 + 504 * s4^3 * (s1)^4 * (s2) * (s3)^2 + 288 * s4^3 * (s1)^2 * (s2)^2 * (s3)^2 + 180 * s4^2 * (s1)^6 * (s2)^2 * (s3)^2 - 96 * s4^3 * (s2)^3 * (s3)^2 + 360 * s4^2 * (s1)^4 * (s2)^3 * (s3)^2 + 144 * s4^2 * (s2)^5 * (s3)^2 + 252 * s4 * (s1)^4 * (s2)^5 * (s3)^2 - 72 * s4 * (s2)^7 * (s3)^2 + 54 * (s1)^2 * (s2)^8 * (s3)^2 + 12 * (s2)^9 * (s3)^2 - 512 * s4^3 * (s1)^3 * (s3)^3 - 24 * s4^2 * (s1)^7 * (s3)^3 - 192 * s4^3 * (s1) * (s2) * (s3)^3 - 432 * s4^2 * (s1)^5 * (s2) * (s3)^3 - 576 * s4^2 * (s1)^3 * (s2)^2 * (s3)^3 - 192 * s4^2 * (s1) * (s2)^3 * (s3)^3 - 240 * s4 * (s1)^5 * (s2)^3 * (s3)^3 - 480 * s4 * (s1)^3 * (s2)^4 * (s3)^3 + 144 * s4 * (s1) * (s2)^5 * (s3)^3 - 112 * (s1)^3 * (s2)^6 * (s3)^3 - 96 * (s1) * (s2)^7 * (s3)^3 + 48 * s4^3 * (s3)^4 + 342 * s4^2 * (s1)^4 * (s3)^4 + 504 * s4^2 * (s1)^2 * (s2) * (s3)^4 + 60 * s4 * (s1)^6 * (s2) * (s3)^4 + 72 * s4^2 * (s2)^2 * (s3)^4 + 540 * s4 * (s1)^4 * (s2)^2 * (s3)^4 + 360 * s4 * (s1)^2 * (s2)^3 * (s3)^4 - 60 * s4 * (s2)^4 * (s3)^4 + 105 * (s1)^4 * (s2)^4 * (s3)^4 + 252 * (s1)^2 * (s2)^5 * (s3)^4 + 42 * (s2)^6 * (s3)^4 - 144 * s4^2 * (s1) * (s3)^5 - 72 * s4 * (s1)^5 * (s3)^5 - 432 * s4 * (s1)^3 * (s2) * (s3)^5 - 144 * s4 * (s1) * (s2)^2 * (s3)^5 - 36 * (s1)^5 * (s2)^2 * (s3)^5 - 240 * (s1)^3 * (s2)^3 * (s3)^5 - 180 * (s1) * (s2)^4 * (s3)^5 + 120 * s4 * (s1)^2 * (s3)^6 + 2 * (s1)^6 * (s3)^6 + 24 * s4 * (s2) * (s3)^6 + 60 * (s1)^4 * (s2) * (s3)^6 + 180 * (s1)^2 * (s2)^2 * (s3)^6 + 40 * (s2)^3 * (s3)^6 - 24 * (s1)^3 * (s3)^7 - 48 * (s1) * (s2) * (s3)^7 + 3 * (s3)^8
    
    return ZZ['x']([n4,-n3,n2,-n1,1])
