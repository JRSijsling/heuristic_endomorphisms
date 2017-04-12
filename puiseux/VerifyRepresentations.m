/***
 *  Verifying representations
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


intrinsic RepresentationCertificate(X::Crv, A::., P::Pt) -> .
{Gives certificate.}

if IsScalar(A) then
    return "Scalar: OK";
else
    // TODO: Add theoretical lower bounds (but for now we can experiment)
    LowerBound := 1;
    return DivisorFromMatrixSplit(X, P, Transpose(A) : LowerBound := LowerBound);
end if;

end intrinsic;


intrinsic VerifyRepresentations(X::Crv, As::SeqEnum, P::Pt) -> List
{Verifies representations and gives certificates.}

return true, [* RepresentationCertificate(X, A, P) : A in As *];

end intrinsic;
