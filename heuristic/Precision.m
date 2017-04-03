/***
 *  Relative number field functionality
 *
 *  Copyright (C) 2016, 2017 Edgar Costa, Jeroen Sijsling
 *                                       (jeroen.sijsling@uni-ulm.de)
 *  See LICENSE.txt for license details.
 */


declare attributes FldCom : epscomp, epsLLL, epsinv;
declare attributes FldRe  : epscomp, epsLLL, epsinv;


intrinsic ComplexFieldExtra(prec::RngIntElt) -> FldCom
{Creates a complex field with some extra needed parameters.}

CC := ComplexField(prec);
RR := RealField(CC);
CC`epscomp := CC ! (10^(-prec + 30));
CC`epsLLL  := CC ! (5^(-prec + 7));
CC`epsinv  := CC ! (2^(-prec + 30));
RR`epscomp := RR ! (10^(-prec + 30));
RR`epsLLL  := RR ! (5^(-prec + 7));
RR`epsinv  := RR ! (2^(-prec + 30));
return CC;

end intrinsic;


intrinsic SetEpscomp(CC::FldCom, epscomp::.)
{Modifies epscomp.}

RR := RealField(CC);
CC`epscomp := CC ! epscomp;
RR`epscomp := RR ! epscomp;

end intrinsic;


intrinsic SetEpsLLL(CC::FldCom, epsLLL::.)
{Modifies epslll.}

RR := RealField(CC);
CC`epsLLL := CC ! epsLLL;
RR`epsLLL := RR ! epsLLL;

end intrinsic;


intrinsic SetEpsinv(CC::FldCom, epsinv::.)
{Modifies epsinv.}

RR := RealField(CC);
CC`epsinv := CC ! epsinv;
RR`epsinv := RR ! epsinv;

end intrinsic;
