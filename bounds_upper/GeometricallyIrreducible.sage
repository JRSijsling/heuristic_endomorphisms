attach("constants.sage")

def IsGeometricallyIrreducible(LPolys, conductor):
    for p in range (2,maxP):                                # loop over primes
      if is_prime(p) and conductor % p <> 0 :               # ensure p is of good reduction
       # print "Testing p =",p
       q = LPolys[p]                                        
       g = q.degree() / 2
       q = twistPolynomial(q, extensionBounds[g])
       if(q.coefficients(sparse=false)[g] % p<>0) :         # check for ordinarity
            if(q.is_irreducible()) :
                return true;
    return false;
