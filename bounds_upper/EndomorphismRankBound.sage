attach("constants.sage");
attach("DiscriminantBound.sage");
attach("NonQM.sage");
attach("GeometricallyIrreducible.sage");

def EndomorphismRankBound ( LPolys, conductor, genus, provenReducible = false, provenQM = false ) :
    if genus > 3 :
        raise ValueError("The genus is too large")
        
    if genus==2 and provenQM :
        return 4;                                                   # quaternion algebra case

    Irreducible = IsGeometricallyIrreducible(LPolys, conductor)
    
    if not Irreducible and genus == 3 :
        raise NotImplementedError("Reducible threefold")
    
    if Irreducible and provenReducible :
        raise ValueError("The Jacobian is geometrically irreducible, but provenReducible is set to True")
        
    if Irreducible :
        Type, DiscBound = DiscriminantBound( LPolys, conductor )
        if Type == "Z" :
            return 1;
        if Type == "Quadratic" :
            return 2;
        if Type == "RM" :
            return genus;
        if Type == "FullCM" :
            return 2*genus;

    raise NotImplementedError()
