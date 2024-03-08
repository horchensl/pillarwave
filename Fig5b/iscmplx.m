function y = iscmplx(x)
%ISCMPLX ISCMPLX(X) returns a 0 if all imaginary
%	 parts of matrix X are zero, and a 1 if
%TSMAT	 any of them is nonzero.
%TLBX	 ISCMPLX is the complement of ISREAL.
%
%	 See also ISREAL.

%	Uilke Stelwagen, June 1991.
%	Copyright (C) 1992, Institute of Applied Physics, TNO-TPD,
%	The Netherlands.

y=any(imag(x));
