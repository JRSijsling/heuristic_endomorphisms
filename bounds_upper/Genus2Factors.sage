def Genus2FactorTwistedLPolys(LPolys) : # the input LPolys need to be in genus 3!
    Genus2LPolys = [0 for k in range(0,maxP)]
    for p in range (2,maxP) :
      if is_prime(p) and LPolys[p] <> 0 :
        q = LPolys[p]
        # print q
        # print "Twisting the L-poly"
        q = twistPolynomial(q, extensionBounds[3])
        # print "Prime", p, "Polynomial", q;

        pieces = q.factor()
        for piece in pieces :
            if piece[0].degree() == 4 :
                Genus2LPolys[p] = piece[0];

    return Genus2LPolys
