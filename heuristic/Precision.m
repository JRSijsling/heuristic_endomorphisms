/***
 *  Relative number field functionality
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


declare attributes FldCom : epscomp, epsLLL, epsinv;


intrinsic ComplexFieldExtra(prec::RngIntElt) -> FldCom
{Creates a complex field with some extra needed parameters.}

CC := ComplexField(prec);
CC`epscomp := CC ! (10^(-prec + 30));
CC`epsLLL  := CC ! (5^(-prec + 7));
CC`epsinv  := CC ! (2^(-prec + 30));
return CC;

end intrinsic;


intrinsic SetEpscomp(CC::FldCom, epscomp::.)
{Modifies epscomp.}

CC`epscomp := CC ! epscomp;

end intrinsic;


intrinsic SetEpsLLL(CC::FldCom, epsLLL::.)
{Modifies epslll.}

CC`epsLLL := CC ! epsLLL;

end intrinsic;


intrinsic SetEpsinv(CC::FldCom, epsinv::.)
{Modifies epsinv.}

CC`epsinv := CC ! epsinv;

end intrinsic;
