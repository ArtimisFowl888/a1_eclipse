;This code is modified from a1.pro
;Only read_hrrd is modified to accomodate new data format:
;ftp://sparc-ftp1.ceda.ac.uk/sparc/hres/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; PURPOSE : get the projection (i.e., array index) of array B onto array A
;           assuming that B is a subset of A
;
;           more general than project.pro and project2.pro
;
; INPUT : a (bigger array), b(smaller array)
;
; OUTPUT : ava (array index)
;
; ALGORITHM :
;
; NOTES : created on 7/25/2004, more general than projection.pro.old3
;
;	  updated on 10/19/04 to add the loose option
;	  updated on 03/09/05 to add the string keyword
;	  updated on 10/23/05 to add the negative keyword
;	  updated on 10/23/05 to add the sloppy keyword
;
;	  updated on 3/25/21 to adjust the reader to fit eclipse2020 format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;assume B is a subset of A
;B .le. A
function projection,big $		;1-D big (relative to small) array
                   ,small $		;1-D small (relative to big) array
                   ,loose=loose $	;a(i) and b(j) are considered the same
					;if |a(i)-b(j)| <= loose, default:1e-38
                   ,string=string $	;for string type
                   ,negative=negative $ ;negative projection of B on A
                   ,sloppy=sloppy $ 	;B is a not necessary a subset of A
                   ,count=count 	;B is a not necessary a subset of A
on_error,2
;----------------------------------------------
if not keyword_set(loose) and not keyword_set(string) then loose=1.e-38
a=big & b=small & na=n_elements(a) & nb=n_elements(b)
if na lt nb and not keyword_set(sloppy) then $
    message,'Usage: ava=projection,a,b (size(a) ge size(b))'
;----------------------------------------------	negative 0
if not keyword_set(negative) then begin
    ava=lonarr(nb)
    ;------------------------------------------	string 0
    if not keyword_set(string) then begin
        ;--------------------------------------	sloppy 0
        if not keyword_set(sloppy) then begin
            for i=0,nb-1 do begin
                a2=where(abs(b(i)-a) le loose,nava)
                if nava eq 0 then message,'B is not a subset of A'
                ava(i)=a2(0)
            endfor
        ;--------------------------------------	sloppy 1
        endif else begin
            j=0L
            for i=0,nb-1 do begin
                a2=where(abs(b(i)-a) le loose,nava)
                if nava eq 0 then begin
                    print,b(i),' in B is not a subset of A'
                endif else begin
                    ava(j)=a2(0)
                    j=j+1L
                endelse
            endfor
            ava=ava(0:j-1)
        endelse
    ;------------------------------------------	string 1
    endif else begin
        ;--------------------------------------	sloppy 0
        if not keyword_set(sloppy) then begin
            for i=0,nb-1 do begin
                a2=where(b(i) eq a,nava)
                if nava eq 0 then message,'B is not a subset of A'
                ava(i)=a2(0)
            endfor
        ;--------------------------------------	sloppy 1
        endif else begin
            j=0L
            for i=0,nb-1 do begin
                a2=where(b(i) eq a,nava)
                if nava eq 0 then begin
                    print,b(i), ' in B is not a subset of A'
                endif else begin
                    ava(j)=a2(0) & j=j+1L
                endelse
            endfor
            ava=ava(0:j-1)
        endelse
    endelse
;----------------------------------------------	negative 1
endif else begin
    if na eq nb then begin
        print,'size(a) eq size(b) and /negative set'
        count=0
        return,-999.
    endif
    ava=lonarr(na) & j=0L
    ;------------------------------------------	string 0
    if not keyword_set(string) then begin
        for i=0,na-1 do begin
            a2=where(abs(b-a(i)) lt loose,nava)
            if nava eq 0 then begin
                if j eq na-nb then begin
                    if not keyword_set(sloppy) then $
                        message,'B is not a subset of A' $
                    else $
                        print,'B is not a subset of A'
                endif
                ava(j)=i & j=j+1L
            endif
        endfor
        ava=ava(0:j-1)
    ;------------------------------------------	string 1
    endif else begin
        for i=0,na-1 do begin
            a2=where(a(i) eq b,nava)
            if nava eq 0 then begin
                if j eq na-nb then begin
                    if not keyword_set(sloppy) then $
                        message,'B is not a subset of A' $
                    else $
                        print,'B is not a subset of A'
                endif
                ava(j)=i & j=j+1L
            endif
        endfor
        ava=ava(0:j-1)
    endelse
endelse
;----------------------------------------------
count=n_elements(ava)
return,ava
end
pro check,Test_variable,ttl=ttl

;+
; NAME:
;  check
; PURPOSE:
;  does a help of the variable, and if it is an array
;  prints out the max and min
; CATEGORY:
;  debug
; CALLING SEQUENCE:
;   check,test_variable
; INPUTS:
;   test_variable = the variable to be inspected
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
; OPTIONAL OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; REQUIRED ROUTINES:
; MODIFICATION HISTORY:
;   mr schoeberl '88
;  idlv2 (lrlait) 900615
;-

on_error,2

cc=!c

v=size(Test_variable)
if v(0) eq 0 then begin
   help,Test_variable
endif else begin
   if not keyword_set(ttl) then begin
       help,Test_variable
       print,' Min =',min(Test_variable),' Max = ',max(Test_variable)
   endif else begin
       print,ttl+' Min =',min(Test_variable),' Max = ',max(Test_variable)
   endelse
endelse

!c=cc

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;purpose:       cal. k (zonal wavenumber)
;                    l (meridional wavenumber)
;               from
;                    phi (horizontal phase propagation direction)
;                    l_h (horizontal wavelength)
;
;input:
;
;output:
;
;algorithm:
;
;history:       09/16/04 created
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro phi2kl $

          ;input
          ,phi $        ;horizontal phase prop. dir., from East         [deg]
          ,l_h $        ;horizontal wavelength                          [km]

          ;output
          ,k $          ;zonal wavenumber                               [SI]
          ,l            ;meridional wavenumber                          [SI]

;kl=2.*!pi*1.e-3/l_h
;phi2=phi/180.*!pi
kl=2.e-3*!pi/l_h
phi2=phi*!dtor
k=kl*cos(phi2)
l=kl*sin(phi2)

return
end

FUNCTION Avg1D, A0, NoData=NoData, MinPts=MinPts
;+
; PURPOSE:
;       Returns average of 1-d array of data.
;
; CALLING SEQUENCE:
;       Avg1D, A, NoData=NoData, MinPts=MinPts
;
; INPUTS:
;       A:              input vector
;       NoData:         missing data value
;       MinPts:         minimum number of valid data in A (default: 1)
;
; OUTPUTS:
;       Returns TOTAL(WHERE(A NE NoData))/N_ELEMENTS(WHERE(A NE NoData))
;       if N_ELEMENTS(...) is greater then MinPts; returns NoData otherwise.
;_
a=reform(a0)
 IF KEYWORD_SET(MinPts) THEN _MinPts=MinPts ELSE _MinPts=1
 IF KEYWORD_SET(NoData) THEN BEGIN
   in = WHERE(A NE NoData, NIn)
   IF (NIn GE _MinPts) THEN Value=TOTAL(A(in))/NIn ELSE Value=NoData
 ENDIF ELSE Value=TOTAL(A)/N_ELEMENTS(A)
RETURN, Value
END; {Avg1D}

;-------------------------------------------------------------
;+
; NAME:
;       DATATYPE
; PURPOSE:
;       Datatype of variable as a string (3 char or spelled out).
; CATEGORY:
; CALLING SEQUENCE:
;       typ = datatype(var, [flag])
; INPUTS:
;       var = variable to examine.         in
;       flag = output format flag (def=0). in
; KEYWORD PARAMETERS:
;       Keywords:
;         /DESCRIPTOR returns a descriptor for the given variable.
;           If the variable is a scalar the value is returned as
;           a string.  If it is an array a description is return
;           just like the HELP command gives.  Ex:
;           datatype(fltarr(2,3,5),/desc) gives
;             FLTARR(2,3,5)  (flag always defaults to 3 for /DESC).
; OUTPUTS:
;       typ = datatype string or number.   out
;          flag=0    flag=1      flag=2    flag=3
;          UND       Undefined   0         UND
;          BYT       Byte        1         BYT
;          INT       Integer     2         INT
;          LON       Long        3         LON
;          FLO       Float       4         FLT
;          DOU       Double      5         DBL
;          COM       Complex     6         COMPLEX
;          STR       String      7         STR
;          STC       Structure   8         STC
;          DCO       DComplex    9         DCOMPLEX
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       Written by R. Sterner, 24 Oct, 1985.
;       RES 29 June, 1988 --- added spelled out TYPE.
;       R. Sterner, 13 Dec 1990 --- Added strings and structures.
;       R. Sterner, 19 Jun, 1991 --- Added format 3.
;       R. Sterner, 18 Mar, 1993 --- Added /DESCRIPTOR.
;       R. Sterner, 1995 Jul 24 --- Added DCOMPLEX for data type 9.
;       Johns Hopkins University Applied Physics Laboratory.
;
; Copyright (C) 1985, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	function datatype,var, flag0, descriptor=desc, help=hlp
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Datatype of variable as a string (3 char or spelled out).'
	  print,' typ = datatype(var, [flag])'
	  print,'   var = variable to examine.         in'
	  print,'   flag = output format flag (def=0). in'
	  print,'   typ = datatype string or number.   out'
	  print,'      flag=0    flag=1      flag=2    flag=3'
	  print,'      UND       Undefined   0         UND'
	  print,'      BYT       Byte        1         BYT'
	  print,'      INT       Integer     2         INT'
	  print,'      LON       Long        3         LON'
	  print,'      FLO       Float       4         FLT'
	  print,'      DOU       Double      5         DBL'
	  print,'      COM       Complex     6         COMPLEX'
	  print,'      STR       String      7         STR'
	  print,'      STC       Structure   8         STC'
	  print,'      DCO       DComplex    9         DCOMPLEX'
	  print,' Keywords:'
	  print,'   /DESCRIPTOR returns a descriptor for the given variable.'
 	  print,'     If the variable is a scalar the value is returned as'
 	  print,'     a string.  If it is an array a description is return'
 	  print,'     just like the HELP command gives.  Ex:'
 	  print,'     datatype(fltarr(2,3,5),/desc) gives'
 	  print,'       FLTARR(2,3,5)  (flag always defaults to 3 for /DESC).'
	  return, -1
	endif 
 
	if n_params(0) lt 2 then flag0 = 0	; Default flag.
	flag = flag0				; Make a copy.
 
	if n_elements(var) eq 0 then begin
	  s = [0,0]
	endif else begin
	  s = size(var)
	endelse
 
	if keyword_set(desc) then flag = 3
 
	if flag eq 2 then typ = s(s(0)+1)
 
	if flag eq 0 then begin
	  case s(s(0)+1) of
   0:	    typ = 'UND'
   1:       typ = 'BYT'
   2:       typ = 'INT'
   4:       typ = 'FLO'
   3:       typ = 'LON'
   5:       typ = 'DOU'
   6:       typ = 'COM'
   7:       typ = 'STR'
   8:       typ = 'STC'
   9:       typ = 'DCO'
else:       print,'Error in datatype'
	  endcase
	endif else if flag eq 1 then begin
	  case s(s(0)+1) of
   0:	    typ = 'Undefined'
   1:       typ = 'Byte'
   2:       typ = 'Integer'
   4:       typ = 'Float'
   3:       typ = 'Long'
   5:       typ = 'Double'
   6:       typ = 'Complex'
   7:       typ = 'String'
   8:       typ = 'Structure'
   9:       typ = 'DComplex'
else:       print,'Error in datatype'
	  endcase
	endif else if flag eq 3 then begin
	  case s(s(0)+1) of
   0:	    typ = 'UND'
   1:       typ = 'BYT'
   2:       typ = 'INT'
   4:       typ = 'FLT'
   3:       typ = 'LON'
   5:       typ = 'DBL'
   6:       typ = 'COMPLEX'
   7:       typ = 'STR'
   8:       typ = 'STC'
   9:       typ = 'DCOMPLEX'
else:       print,'Error in datatype'
	  endcase
	endif
 
	if ~ keyword_set(desc) then begin
	  return, typ					; Return data type.
	endif else begin
	  if s(0) eq 0 then return,strtrim(var,2)	; Return scalar desc.
	  aa = typ+'ARR('
          for i = 1, s(0) do begin                      
            aa = aa + strtrim(s(i),2)                 
            if i lt s(0) then aa = aa + ','          
            endfor                                     
          aa = aa+')'                                   
	  return, aa
	endelse
 
	END
;---------------------------------------------------------------------------
PRO  resetplt, all=all,x=x,y=y,z=z,p=p, invert=invert
;+
; NAME:			RESETPLT
;
; PURPOSE:		This procedure will reset all or a selection 
;			of the system plot variables to their initial values
;
; CATEGORY:		Plot Utility
;
; CALLING SEQUENCE:	resetplt,/all		;clear the !p, !x, !y, !z 
;			resetplt,/x,/z		;clear the !x and !z variables 
;			resetplt,/x		;only clear the !x variable
;			resetplt,/x,/invert	;clear all except the !x 
;
; INPUTS:		
;	KEYWORDS:
;		x,y,z,p	= clear the appropriate variable
;		all	= clear all, this is equivalent to /x,/y,/z,/p
;		invert	= invert the logic. Clear all unselected variables.
;			  Therefore "clearplt,/all,/invert" does nothing.
;
; OUTPUTS:	none
;
; COMMON BLOCKS:
;	none.
; SIDE EFFECTS:
;		The sytem plot variables are changed.
;	
; MODIFICATION HISTORY:
;	Written by: Trevor Harris, Physics Dept., University of Adelaide,
;		July, 1990.
;
;-

   resetx = 0b
   resety = 0b
   resetz = 0b
   resetp = 0b
   if (keyword_set(x)) then resetx = 1b
   if (keyword_set(y)) then resety = 1b
   if (keyword_set(z)) then resetz = 1b
   resetp = not (resetx or resety or resetz) 
   if (keyword_set(p)) then resetp = 1b
   if (keyword_set(all)) then begin
      resetx = 1b
      resety = 1b
      resetz = 1b
      resetp = 1b
   endif

   if (keyword_set(invert)) then begin
      resetx = not resetx
      resety = not resety
      resetz = not resetz
      resetp = not resetp
   endif

   if (resetx) then begin
      !x.thick=0.0
      !x.charsize=0.0
      !x.ticks=0
      !x.tickv=0
      !x.tickname=''
      !x.title=' '
      !x.range=0
      !x.ticklen=0
      !x.style=0
      !x.margin = [10,3]
      !x.tickformat=''
   endif
   if (resety) then begin
      !y.thick=0.0
      !y.charsize=0.0
      !y.ticks=0
      !y.tickv=0
      !y.tickname=''
      !y.title=' '
      !y.range=0
      !y.ticklen=0
      !y.style=0
      !y.margin = [4,2]
      !y.tickformat=''
   endif
   if (resetz) then begin
      !z.thick=0.0
      !z.charsize=0.0
      !z.ticks=0
      !z.tickv=0
      !z.tickname=''
      !z.title=' '
      !z.range=0
      !z.ticklen=0
      !z.style=0
      !z.margin = [0,0]
      !z.tickformat=''
   endif
   if (resetp) then begin
      !p.title=' '
      !p.subtitle=' '
      !p.ticklen=0.02
      !p.charsize=1.0
      !p.charthick=1.0
      !p.thick=1.0
      !p.linestyle=0
      !p.region = [0,0,0,0]
      !p.position = [0,0,0,0]
      !p.psym = 0
      !p.nsum = 0
   endif

   return
END 
;-------------------------------------------------------------
;+
; NAME:
;       YMD2JD
; PURPOSE:
;       From Year, Month, and Day compute Julian Day number.
; CATEGORY:
; CALLING SEQUENCE:
;       jd = ymd2jd(y,m,d)
; INPUTS:
;       y = Year (like 1987).                    in
;       m = month (like 7 for July).             in
;       d = month day (like 23).                 in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       jd = Julian Day number (like 2447000).   out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner,  23 June, 1985 --- converted from FORTRAN.
;       Johns Hopkins University Applied Physics Laboratory.
;       RES 18 Sep, 1989 --- converted to SUN
;
; Copyright (C) 1985, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	function ymd2jd, iy, im, id, help=hlp
 
	if (n_params(0) LT 3) or keyword_set(hlp) then begin
	  print,' From Year, Month, and Day compute Julian Day number.'
	  print,' jd = ymd2jd(y,m,d)'
	  print,'   y = Year (like 1987).                    in'
	  print,'   m = month (like 7 for July).             in'
	  print,'   d = month day (like 23).                 in'
	  print,'   jd = Julian Day number (like 2447000).   out'
	  return, -1
	endif
 
	y = long(iy)
	m = long(im)
	d = long(id)
	jd = 367*y-7*(y+(m+9)/12)/4-3*((y+(m-9)/7)/100+1)/4 $
             +275*m/9+d+1721029
 
	return, jd
 
	end
;-------------------------------------------------------------
;+
; NAME:
;       JS2JD
; PURPOSE:
;       Convert from Julian Seconds to Julian Day Number.
; CATEGORY:
; CALLING SEQUENCE:
;       jd = js2jd(js)
; INPUTS:
;       js = Equivalent Julian Second.    in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       jd = Julian Day Number.           out
; COMMON BLOCKS:
; NOTES:
;       Notes: JS are seconds after 2000 Jan 1 0:00.
; MODIFICATION HISTORY:
;       R. Sterner, 1994 Oct 13
;
; Copyright (C) 1994, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	function js2jd, js, help=hlp
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Convert from Julian Seconds to Julian Day Number.'
	  print,' jd = js2jd(js)'
 	  print,'   js = Equivalent Julian Second.    in'
	  print,'   jd = Julian Day Number.           out'
	  print,' Notes: JS are seconds after 2000 Jan 1 0:00.'
	  return,''
	endif
 
	jd2000 = ymd2jd(2000,1,1) - 0.5d0	; JD at 2000 Jan 1 0:00.
	jd = jd2000 + js/86400d0		; JD for JS.
 
	return, jd
	end
;-------------------------------------------------------------
;+
; NAME:
;       JS2YMDS
; PURPOSE:
;       Convert from "Julian Second" to year, month, day, second.
; CATEGORY:
; CALLING SEQUENCE:
;       js2ymds, js, y, m, d, s
; INPUTS:
;       js = "Julian Second".               in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       y,m,d = year, month, day numbers.   out
;       s = second into day.                out
; COMMON BLOCKS:
; NOTES:
;       Notes: Julian seconds (not an official unit) serve the
;         same purpose as Julian Days, interval computations.
;         The zero point is 0:00 1 Jan 2000, so js < 0 before then.
;         Julian Seconds are double precision and have a precision
;         better than 1 millisecond over a span of +/- 1000 years.
;         A precision warning may point to a call to dt_tm_fromjs.
;
;       See also ymds2js, dt_tm_tojs, dt_tm_fromjs, jscheck.
; MODIFICATION HISTORY:
;       R. Sterner, 2 Sep, 1992
;       R. Sterner, 13 Dec, 1992 --- added data type check.
;
; Copyright (C) 1992, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------

        pro js2ymds, js, y, m, d, s, help=hlp

        if (n_params(0) lt 4) or keyword_set(hlp) then begin
          print,' Convert from "Julian Second" to year, month, day, second.'
          print,' js2ymds, js, y, m, d, s'
          print,'   js = "Julian Second".               in'
          print,'   y,m,d = year, month, day numbers.   out'
          print,'   s = second into day.                out'
          print,' Notes: Julian seconds (not an official unit) serve the'
          print,'   same purpose as Julian Days, interval computations.'
          print,'   The zero point is 0:00 1 Jan 2000, so js < 0 before then.'
          print,'   Julian Seconds are double precision and have a precision'
          print,'   better than 1 millisecond over a span of +/- 1000 years.'
          print,'   A precision warning may point to a call to dt_tm_fromjs.'
          print,' '
          print,' See also ymds2js, dt_tm_tojs, dt_tm_fromjs, jscheck.'
          return
        endif

        sz = size(js)
        if sz(sz(0)+1) ne 5 then begin
          print,' Warning in js2ymds: Julian Seconds should be passed in'
          print,'   as double precision.  Precision degraded.'
        endif
        days = floor(js/86400)
        s = js - 86400D0*days
        jd2ymd, days+2451545, y, m, d

        return
        end
;-------------------------------------------------------------
;+
; NAME:
;       JD2YMD
; PURPOSE:
;       Find year, month, day from julian day number.
; CATEGORY:
; CALLING SEQUENCE:
;       jd2ymd, jd, y, m, d
; INPUTS:
;       jd = Julian day number (like 2447000).     in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       y = year (like 1987).                      out
;       m = month number (like 7).                 out
;       d = day of month (like 23).                out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner.  21 Aug, 1986.
;       Johns Hopkins Applied Physics Lab.
;       RES 18 Sep, 1989 --- converted to SUN
;       R. Sterner, 30 Apr, 1993 --- cleaned up and allowed arrays.
;
; Copyright (C) 1986, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------

        pro jd2ymd, jdd, y, m, d, help=hlp

        if (n_params(0) lt 4) or keyword_set(hlp) then begin
          print,' Find year, month, day from julian day number.'
          print,' jd2ymd, jd, y, m, d'
          print,'   jd = Julian day number (like 2447000).     in'
          print,'   y = year (like 1987).                      out'
          print,'   m = month number (like 7).                 out'
          print,'   d = day of month (like 23).                out'
          return
        endif

        jd = long(jdd)                          ; Force long.
        y = fix((jd - 1721029)/365.25)          ; Estimated year.
        jd0 = ymd2jd(y, 1, 0)                   ; JD for day 0.
        days = jd - jd0                         ; Day of year.
        w = where(days le 0, cnt)               ; Find where year is wrong.
        if cnt gt 0 then begin
          y(w) = y(w) - 1                       ; Year was off by 1.
          jd0(w) = ymd2jd( y(w), 1, 0)          ; New JD for day 0.
          days(w) = jd(w) - jd0(w)              ; New day of year.
        endif

        ;---  Correct for leap-years.  -----
        ly = (((y mod 4) eq 0) and ((y mod 100) ne 0)) $
            or ((y mod 400) eq 0)

        ;---  Days before start of each month.  -----
        ydays = [0,0,31,59,90,120,151,181,212,243,273,304,334,366]
        off   = [0,0, 0, 1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1]

        ;----------------  Find which month.  --------------------------
        ;     Algorithm: ydays has cumulative # days up to start of each month
        ;     (13 elements so month number may be used as an index).
        ;     This number needs 1 added for Mar to Dec if it is a leap year.
        ;     This is done by adding the offset, off, times the leap year flag.
        ;     The larger the day of year (days) the fewer elements of this
        ;     corrected ydays array will be greater than days.  The
        ;     entries in the corrected ydays gt days are located by where and
        ;     counted by n_elements.  Ydays has 13 elements so subtract result
        ;     from 13 to get month number.
        ;---------------------------------------------------------------
        njd = n_elements(jd)    ; # of JDs to convert.
        m = intarr(njd)         ; Set up storage for months.
        d = intarr(njd)         ; Set up storage for day of month.
        for i = 0, njd-1 do begin       ; Loop through each JD.
          ydays2 = ydays+ly(i)*off      ; Correct cumulative days for year.
          dy = days(i)                  ; Days into year for i'th JD.
          mn = 13-n_elements(where(dy le ydays2))  ; i'th month number.
          m(i) = mn                     ; Store month.
          d(i) = fix(dy - ydays2(mn))   ; Find and store i'th day of month.
        endfor

        ;---------  Make sure scalars are returned as scalars  -------
        if n_elements(m) eq 1 then begin
          m = m(0)
          d = d(0)
        endif

        return
        end
;-------------------------------------------------------------
;+
; NAME:
;       YMDS2JS
; PURPOSE:
;       Convert to year, month, day, second to "Julian Second".
; CATEGORY:
; CALLING SEQUENCE:
;       js = ymds2js(y,m,d,s)
; INPUTS:
;       y,m,d = year, month, day numbers.   in
;       s = second into day.                in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       js = "Julian Second".               out
; COMMON BLOCKS:
; NOTES:
;       Notes: Julian seconds (not an official unit) serve the
;         same purpose as Julian Days, interval computations.
;         The zero point is 0:00 1 Jan 2000, so js < 0 before then.
;         Julian Seconds are double precision and have a precision
;         better than 1 millisecond over a span of +/- 1000 years.
;
;       See also js2ymds, dt_tm_fromjs, dt_tm_tojs, jscheck.
; MODIFICATION HISTORY:
;       R. Sterner, 2 Sep, 1992
;
; Copyright (C) 1992, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------

        function ymds2js, y, m, d, s, help=hlp

        if (n_params(0) lt 4) or keyword_set(hlp) then begin
          print,' Convert to year, month, day, second to "Julian Second".'
          print,' js = ymds2js(y,m,d,s)'
          print,'   y,m,d = year, month, day numbers.   in'
          print,'   s = second into day.                in'
          print,'   js = "Julian Second".               out'
          print,' Notes: Julian seconds (not an official unit) serve the'
          print,'   same purpose as Julian Days, interval computations.'
          print,'   The zero point is 0:00 1 Jan 2000, so js < 0 before then.'
          print,'   Julian Seconds are double precision and have a precision'
          print,'   better than 1 millisecond over a span of +/- 1000 years.'
          print,' '
          print,' See also js2ymds, dt_tm_fromjs, dt_tm_tojs, jscheck.'
          return, -1
        endif

        return, s + (ymd2jd(y,m,d)-2451545)*86400d0

        end
;********************************
PRO arc_circ,r,theta_s,df
; Procedure to plot arc of circle with radius r,
; starting angle theta_s over and arc length df
   x = fltarr(100)
   y = x
   x = replicate(r,100)
   y = findgen(100)*df/100+theta_s
   oplot,x,y,/polar
   xo = [0.,x*cos(y)]
   yo = [0.,x*sin(y)]
   polyfill,xo,yo,color = 200
   return
END
;-------------------------------------------------------------
;+
; NAME:
;       INTERPX
; PURPOSE:
;       Interpolate data with possible gaps and missing (bad) values.
; CATEGORY:
; CALLING SEQUENCE:
;       yy = interp1(x,y,xx)
; INPUTS:
;       x,y = input points.            in
;         x is assumed monotonically increasing.
;       xx = desired x values.         in
;         xx need not be monotonically increasing.
; KEYWORD PARAMETERS:
;       Keywords:
;         BAD=b  Values GREATER THAN b are considered missing.
;         GAP=g  Gaps in x greater than or equal to this
;           are broken by setting the output curve points
;           in the gaps to a flag value of 32000 or BAD if given.
;         /FIXBAD  means interpolate across bad data between
;           closest good points on each side.  Otherwise the
;           bad points are flagged with the value specified
;           for BAD.
; OUTPUTS:
;       yy = interpolated y values.    out
; COMMON BLOCKS:
; NOTES:
;       Notes: Flagged values may be used to break a plotted
;          curve using MAX_VALUE in the PLOT or OPLOT command:
;          plot,x,y,max_value=999
;          SLOW for more than a few thousand points.
; MODIFICATION HISTORY:
;       R. Sterner, 12 Aug, 1993
;
; Copyright (C) 1993, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	function interpx, x,y,xx, help=hlp, bad=bad, $
	  gap=gap, fixbad=xbad
 
	if (n_params(0) lt 3) or keyword_set(hlp) then begin
	  print,' Interpolate data with possible gaps and missing (bad) values.'
	  print,' yy = interp1(x,y,xx)'
	  print,'   x,y = input points.            in'
	  print,'     x is assumed monotonically increasing.'
	  print,'   xx = desired x values.         in'
	  print,'     xx need not be monotonically increasing.'
	  print,'   yy = interpolated y values.    out'
	  print,' Keywords:'
	  print,'   BAD=b  Values GREATER THAN b are considered missing.'
	  print,'   GAP=g  Gaps in x greater than or equal to this'
	  print,'     are broken by setting the output curve points'
	  print,'     in the gaps to a flag value of 32000 or BAD if given.'
	  print,'   /FIXBAD  means interpolate across bad data between'
	  print,'     closest good points on each side.  Otherwise the'
	  print,'     bad points are flagged with the value specified'
	  print,'     for BAD.'
	  print,' Notes: Flagged values may be used to break a plotted'
	  print,'    curve using MAX_VALUE in the PLOT or OPLOT command:'
	  print,'    plot,x,y,max_value=999'
	  print,'    SLOW for more than a few thousand points.'
	  return,-1
	endif
 
	;--------  Find any data gaps  ------------------------------
	;  Data gaps are where the X coordinate jumps by more than a
	;  specified amount.  Save the index of the point just before
	;  the gap in glo and the index of the point just above the
	;  gap in ghi (= next point).
	;------------------------------------------------------------
	cntg = 0				; Assume no gaps.
	if n_elements(gap) ne 0 then begin	; Gap value given.
	  glo = where(x(1:*)-x gt gap, cntg)	; Index of 1st pt below gap.
	  ghi = glo + 1				; Index of 1st pt above gap.
	endif
 
	;--------  Find any bad points  -----------------------------------
	;  It is assumed that bad points have been tagged with a flag value
	;  before calling this routine.  Any points with a value less than
	;  or equal to the flag value are assumed to be good.  The indices
	;  of the good points on either side are saved in lo and hi.
	;-----------------------------------------------------------------
	cntb = 0				; Assume no bad points.
	if n_elements(bad) ne 0 then begin	; Bad tag value given.
	  w = where(y le bad, cnt)		; Indices of good values.
	  if cnt gt 0 then begin		; Found some good points.
	    wb = where(w(1:*)-w gt 1, cntb)	; Look for index jumps > 1.
	  endif
	endif
	if cntb gt 0 then begin	; Found some bad point groups.
	  lo = w(wb)		; Index of ignore window bottom.
	  hi = w(wb+1)		; Index of ignore window top.
	endif
 
	;--------  Setup output y array and flag value  -------
	yy = make_array(n_elements(xx),type=datatype(y,2))
	flag = 32000.
	if n_elements(bad) ne 0 then flag = bad
 
	;-----  Handle bad points  -------------------------------
	;  Flag which points in the output arrays to ignore due to bad pts.
	;---------------------------------------------------------
	if ~ keyword_set(xbad) then begin
	  if cntb gt 0 then begin
	    for i = 0, n_elements(lo)-1 do begin  ; Loop thru bad point gaps.
	      w = where((xx gt x(lo(i))) and (xx lt x(hi(i))), cnt) ; In gap.
	      if cnt gt 0 then yy(w) = flag	  ; Flag those in bad pt gap.
	    endfor
	  endif
	endif
 
	;-----  Handle data gaps  ---------------------------------
	;  Flag which points in the output arrays to ignore due to data gaps.
	;----------------------------------------------------------
	if cntg gt 0 then begin
	  for i = 0, n_elements(glo)-1 do begin	; Loop thru gaps.
	    w = where((xx gt x(glo(i))) and (xx lt x(ghi(i))), cnt)
	    if cnt gt 0 then yy(w) = flag		; Flag those in gap.
	  endfor
	endif
 
	;-----  Linearly interpolate into x,y at xx to get yy --------
	;  Drop any bad points from input curve before interpolating.
	;  Then find input and output curve sizes.
	;  Finally interpolate needed points.
	;-------------------------------------------------------------
	cnt = 0
	if n_elements(bad) ne 0 then begin
	  w = where(y lt bad, cnt)			; Drop bad points
	endif
	if cnt eq 0 then w = lindgen(n_elements(x))	; before interpolation.
	x2 = x(w)
	y2 = y(w)
 
	lstxx = n_elements(xx) - 1		; Number of output points.
	lstx = n_elements(x2) - 1		; Number of good input points.
 
	for i = 0L, lstxx do begin		; Loop through output points.
	  if yy(i) eq flag then goto, skip	; Ignore flagged points.
	  j = (where(x2 ge xx(i)))(0)
	  case 1 of
j eq -1:    yy(i) = y2(lstx)			; After last input x.
j eq  0:    yy(i) = y2(0)			; Before first input x.
   else:    begin				; Must interpolate.
	      m = (y2(j) - y2(j-1))/(x2(j) - x2(j-1))	; Slope.
	      yy(i) = m*(xx(i) - x2(j-1)) + y2(j-1)	; y = m*x + b.
	    end
	  endcase
skip:
	endfor
 
	return, yy
	end
;+
; NAME:
;   MPFIT
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://astrog.physics.wisc.edu/~craigm/idl.html
;
; PURPOSE:
;   Perform Levenberg-Marquardt least-squares minimization (MINPACK-1)
;
; MAJOR TOPICS:
;   Curve and Surface Fitting
;
; CALLING SEQUENCE:
;   parms = MPFIT(MYFUNCT, start_parms, FUNCTARGS=fcnargs, NFEV=nfev,
;                 MAXITER=maxiter, ERRMSG=errmsg, NPRINT=nprint, 
;                 FTOL=ftol, XTOL=xtol, GTOL=gtol, 
;                 STATUS=status, ITERPROC=iterproc, ITERARGS=iterargs,
;                 COVAR=covar, PERROR=perror,
;                 PARINFO=parinfo)
;
; DESCRIPTION:
;
;  MPFIT uses the Levenberg-Marquardt technique to solve the
;  least-squares problem.  In its typical use, MPFIT will be used to
;  fit a user-supplied function (the "model") to user-supplied data
;  points (the "data") by adjusting a set of parameters.
;
;  For example, a researcher may think that a set of observed data
;  points is best modelled with a Gaussian curve.  A Gaussian curve is
;  parameterized by its mean, standard deviation and normalization.
;  MPFIT will, within certain constraints, find the set of parameters
;  which best fits the data.  The fit is "best" in the least-squares
;  sense; that is, the sum of the weighted squared differences between
;  the model and data is minimumized.
;
;  The Levenberg-Marquardt technique is a particular strategy for
;  iteratively searching for the best fit.  This particular
;  implementation is drawn from MINPACK-1 (see NETLIB), and seems to
;  be more robust than routines provided with IDL.  This version
;  allows upper and lower bounding constraints to be placed on each
;  parameter, or the parameter can be held fixed.
;
;  The IDL user-supplied function should return an array of weighted
;  deviations between model and data.  In a typical scientific problem
;  the residuals should be weighted so that each deviate has a
;  gaussian sigma of 1.0.  If X represents values of the independent
;  variable, Y represents a measurement for each value of X, and ERR
;  represents the error in the measurements, then the deviates could
;  be calculated as follows:
;
;    DEVIATES = (Y - F(X)) / ERR
;
;  where F is the analytical function representing the model.  The
;  convenience functions MPFITFUN and MPFITEXPR calculate the deviates
;  for you.  If ERR are the 1-sigma uncertainties in Y, then
;
;    TOTAL( DEVIATES^2 ) 
;
;  will be the total chi-squared value.  MPFIT will minimize the
;  chi-square value.  The values of X, Y and ERR are passed through
;  MPFIT to the user-supplied function via the FUNCTARGS keyword.
;
;  MPFIT does not perform more general optimization tasks.  It is
;  customized, based on MINPACK-1, to the least-squares minimization
;  problem.
;
;  In the search for the best-fit solution, MPFIT calculates
;  derivatives numerically via a finite difference approximation.  The
;  user-supplied function need not calculate the derivatives
;  explicitly.
;   
; INPUTS:
;   MYFUNCT - a string variable containing the name of the function to
;             be minimized.  The function should return the weighted
;             deviations between the model and the data.  It should be
;             declared in the following way (or in some equivalent):
;
;             FUNCTION MYFUNCT, p, X=x, Y=y, ERR=err
;              ; Parameter values are passed in "p"
;              ; Calculation of deviations occurs here
;              model = F(x)
;              return, (y-model)/err
;             END
;
;             The keyword parameters X, Y, and ERR in the example 
;             above are suggestive but not required.  Any parameters
;             can be passed to MYFUNCT by using the FUNCTARGS keyword
;             to MPFIT.  Use MPFITFUN and MPFITEXPR if you need ideas
;             on how to do that.
;
;             In general there are no restrictions on the number of
;             dimensions in X, Y or ERR.  However the deviates *must*
;             be returned in a one-dimensional array, and must have
;             the same type (float or double) as the input arrays.
;
;   START_PARAMS - An array of starting values for each of the
;                  parameters of the model.  The number of parameters
;                  should be fewer than the number of measurements.
;                  Also, the parameters should have the same data type
;                  as the measurements (double is preferred).
;
;                  This parameter is optional if the PARINFO keyword
;                  is used.  See below.  The PARINFO keyword provides
;                  a mechanism to fix or constrain individual
;                  parameters.  If both START_PARAMS and PARINFO are
;                  passed, then the starting *value* is taken from
;                  START_PARAMS, but the *constraints* are taken from
;                  PARINFO.
; 
; INPUT KEYWORD PARAMETERS:
;
;   FUNCTARGS - A structure which contains the parameters to be passed
;               to the user-supplied function specified by MYFUNCT via
;               the _EXTRA mechanism.  This is the way you can pass
;               data to your user-supplied function without using
;               common blocks.
;
;               Consider the following example:
;                if FUNCTARGS = { XVAL:[1.D,2,3], YVAL:[1.D,4,9],
;                                 ERRVAL:[1.D,1,1] }
;                then the user supplied function should be declared
;                like this:
;                PRO MYFUNCT, P, XVAL=x, YVAL=y, ERRVAL=err
;
;               By default, no extra parameters are passed to the
;               user-supplied function.
;
;   MAXITER - The maximum number of iterations to perform.  If the
;             number is exceeded, then the STATUS value is set to 5
;             and MPFIT returns.
;             Default: 200 iterations
;
;   FTOL - a nonnegative input variable. Termination occurs when both
;          the actual and predicted relative reductions in the sum of
;          squares are at most FTOL.  Therefore, FTOL measures the
;          relative error desired in the sum of squares.
;          Default: 1D-10
;
;   XTOL - a nonnegative input variable. Termination occurs when the
;          relative error between two consecutive iterates is at most
;          XTOL. therefore, XTOL measures the relative error desired
;          in the approximate solution.
;          Default: 1D-10
;
;   GTOL - a nonnegative input variable. Termination occurs when the
;          cosine of the angle between fvec and any column of the
;          jacobian is at most GTOL in absolute value. Therefore, GTOL
;          measures the orthogonality desired between the function
;          vector and the columns of the jacobian.
;          Default: 1D-10
;
;   ITERPROC - The name of a procedure to be called upon each NPRINT
;              iteration of the MPFIT routine.  It should be declared
;              in the following way:
;
;              PRO ITERPROC, MYFUNCT, p, iter, FUNCTARGS=fcnargs, $
;                PARINFO=parinfo, QUIET=quiet, ...
;                ; perform custom iteration update
;              END
;         
;              Where MYFUNCT is the user-supplied function to be
;              minimized, P is the current set of model parameters,
;              ITER is the iteration number, and FUNCTARGS are the
;              arguments to be passed to MYFUNCT.  QUIET is set when
;              no textual output should be printed.  See below for
;              documentation of PARINFO.
;
;              In implementation, ITERPROC, can perform updates to the
;              terminal or graphical user interface, to provide
;              feedback while the fit proceeds.  If the fit is to be
;              stopped for any reason, then ITERPROC should set the
;              system variable !ERR to a negative value.  In
;              principle, ITERPROC should probably not modify the
;              parameter values, because it may interfere with the
;              algorithm's stability.  In practice it is allowed.
;
;              Default: an internal routine is used to print the
;                       parameter values.
;
;   NPRINT - The frequency with which ITERPROC is called.  A value of
;            1 indicates that ITERPROC is called with every iteration,
;            while 2 indicates every other iteration, etc.  
;            Default value: 1
;
;   ITERARGS - The keyword arguments to be passed to ITERPROC via the
;              _EXTRA mechanism.  This should be a structure, and is
;              similar in operation to FUNCTARGS.
;              Default: no arguments are passed.
;
;   PARINFO - Provides a mechanism for more sophisticated constraints
;             to be placed on parameter values.  When PARINFO is not
;             passed, then it is assumed that all parameters are free
;             and unconstrained.  No values in PARINFO are modified
;             during the call to MPFIT.
;
;             PARINFO should be an array of structures, one for each
;             parameter.  Each parameter is associated with one
;             element of the array, in numerical order.  The structure
;             should have at least the following entries:
;
;               - VALUE - the starting parameter value (but see
;                         START_PARAMS above).
;
;               - FIXED - a boolean value, whether the parameter is to 
;                         be held fixed or not.  Fixed parameters are
;                         not varied by MPFIT, but are passed on to 
;                         MYFUNCT for evaluation.
;
;               - LIMITED - a two-element boolean array.  If the
;                 first/second element is set, then the parameter is
;                 bounded on the lower/upper side.  A parameter can be
;                 bounded on both sides.
;
;               - LIMITS - a two-element float or double array.  Gives
;                 the parameter limits on the lower and upper sides,
;                 respectively.  Zero, one or two of these values can
;                 be set, depending on the value of LIMITED.
;
;               - STEP - the step size to be used in calculating the
;                 numerical derivatives.  If set to zero, then the
;                 step size is computed automatically.
; 
;             Other tag values can also be given in the structure, but
;             they are ignored.
;
;             Example:
;             parinfo = replicate({value:0.D, fixed:0, $
;                           limited:[0,0], limits:[0.D,0], step:0.D}, 5)
;             parinfo(0).fixed = 1
;             parinfo(1).limited(4) = 1
;             parinfo(1).limits(4)  = 50.D
;             parinfo(*).value = [5.7D, 2.2, 500., 1.5, 2000.]
;
;             A total of 5 parameters, with starting values of 5.7,
;             2.2, 500, 1.5, and 2000 are given.  The first parameter
;             is fixed at a value of 5.7, and the last parameter is
;             constrained to be above 50.
;
;             Default value:  all parameters are free and unconstrained.
;
;   QUIET - set this keyword when no textual output should be printed
;           by MPFIT
;
;   NOCOVAR - set this keyword to prevent the calculation of the
;             covariance matrix before returning (see COVAR)
;
; RETURNS:
;
;   Returns the array of best-fit parameters.
;
; OUTPUT KEYWORD PARAMETERS:
;
;   NFEV - the number of MYFUNCT function evaluations performed.
;
;   ERRMSG - a string error or warning message is returned.
;
;   BESTNORM - the value of the summed squared residuals for the
;              returned parameter values.
;
;   PERROR - The formal 1-sigma errors in each parameter.  If a
;            parameter is held fixed, or if it touches a boundary,
;            then the error is reported as zero.
;
;   COVAR - the covariance matrix for the set of parameters returned
;           by MPFIT.  The matrix is NxN where N is the number of
;           parameters.  The square root of the diagonal elements
;           gives the formal 1-sigma statistical errors on the
;           parameters IF errors were treated "properly" in MYFUNC.
;           Parameter errors are also returned in PERROR.
;
;           To compute the correlation matrix, PCOR, use this:
;           IDL> PCOR = COV * 0
;           IDL> FOR i = 0, n-1 DO FOR j = 0, n-1 DO $
;                PCOR(i,j) = COV(i,j)/sqrt(COV(i,i)*COV(j,j))
;
;           If NOCOVAR is set or MPFIT terminated abnormally, then
;           COVAR is set to a scalar with value !VALUES.D_NAN.
;
;   STATUS - an integer status code is returned.  It can have one of
;            the following values:
;
;	   0  improper input parameters.
;         
;	   1  both actual and predicted relative reductions
;	      in the sum of squares are at most FTOL.
;         
;	   2  relative error between two consecutive iterates
;	      is at most XTOL
;         
;	   3  conditions for STATUS = 1 and STATUS = 2 both hold.
;         
;	   4  the cosine of the angle between fvec and any
;	      column of the jacobian is at most GTOL in
;	      absolute value.
;         
;	   5  the maximum number of iterations has been reached
;         
;	   6  FTOL is too small. no further reduction in
;	      the sum of squares is possible.
;         
;	   7  XTOL is too small. no further improvement in
;	      the approximate solution x is possible.
;         
;	   8  GTOL is too small. fvec is orthogonal to the
;	      columns of the jacobian to machine precision.
;
; EXAMPLE:
;
;   p0 = [5.7D, 2.2, 500., 1.5, 2000.]
;   fa = {X:x, Y:y, ERR:err}
;   p = mpfit('MYFUNCT', p0, functargs=fa)
;
;   Minimizes sum of squares of MYFUNCT.  MYFUNCT is called with the X,
;   Y, and ERR keyword parameters that are given by FUNCTARGS.  The
;   resulting parameter values are returned in p.
;
; REFERENCES:
;
;   MINPACK-1, Jorge More', available from netlib (www.netlib.org).
;   "Optimization Software Guide," Jorge More' and Stephen Wright, 
;     SIAM, *Frontiers in Applied Mathematics*, Number 14.
;
; MODIFICATION HISTORY:
;   Translated from MINPACK-1 in FORTRAN, Apr-Jul 1998, CM
;   Fixed bug in parameter limits (x vs xnew), 04 Aug 1998, CM
;   Added PERROR keyword, 04 Aug 1998, CM
;
;-

FORWARD_FUNCTION mpfit_fdjac2, mpfit_enorm, mpfit_lmpar, mpfit_covar, mpfit

;  Things to do:
;    * optional derivative in user-supplied function

;  Documentation below is mostly original to the MINPACK-1 FORTRAN
;  routines.  Some documentation has been modified to refer to IDL
;  terminology, but not all.

;     **********
;
;     subroutine fdjac2
;
;     this subroutine computes a forward-difference approximation
;     to the m by n jacobian matrix associated with a specified
;     problem of m functions in n variables.
;
;     the subroutine statement is
;
;       PRO        FDJAC2,fcn, x, fvec, step, ulimited, ulimit, 
;                         iflag=iflag, epsfcn=epsfcn, nfev=nfev,
;                         fcnargs=fcnargs, xall=xall, ifree=ifree
;
;     where
;
;       fcn is a string containing the name of a user-supplied
;         function which calculates the m functions.  fcn should
;         be declared in the following way
;
;         FUNCTION fcn, x, ...
;           ; calculations of y values
;           RETURN, y
;         END
;         
;         Extra keyword arguments can be passed via the fcnargs
;         keyword.  The value of !err should not be changed by fcn
;         unless the user wants to terminate execution of fdjac2.  In
;         this case set !err to a negative integer.
;
;       fcnargs is an optional structure containing keyword arguments
;         to be passed to the function fcn, following the _EXTRA
;         convention.  If no extra keyword arguments are to be passed,
;         then fcnargs should remain undefined.
;
;	x is an input array
;
;       fvec is an input array which must contain the functions
;         evaluated at x.
;
;	fjac is an output m by n array which contains the
;	  approximation to the jacobian matrix evaluated at x.
;         m is the number of elements in fvec, n is then number of
;         elements in x.
;
;	iflag is an integer variable which can be used to terminate
;	  the execution of fdjac2. see description of fcn.
;
;	epsfcn is an input variable used in determining a suitable
;	  step length for the forward-difference approximation. this
;	  approximation assumes that the relative errors in the
;	  functions are of the order of epsfcn. if epsfcn is less
;	  than the machine precision, it is assumed that the relative
;	  errors in the functions are of the order of the machine
;	  precision.
;
;	wa is a work array of length m.
;
;     subprograms called
;
;	user-supplied ...... fcn
;
;	minpack-supplied ... dpmpar
;
;	fortran-supplied ... dabs,dmax1,dsqrt
;
;     argonne national laboratory. minpack project. march 1980.
;     burton s. garbow, kenneth e. hillstrom, jorge j. more
;

function mpfit_fdjac2, fcn, x, fvec, step, ulimited, ulimit, $
                 iflag=iflag, epsfcn=epsfcn, nfev=nfev, $
                 FUNCTARGS=fcnargs, xall=xall, ifree=ifree
  sz = size(x)
  mch = machar(double=(sz(sz(0)+1) EQ 5))
  MACHEP = mch.eps
  DWARF = mch.xmin
  if n_elements(epsfcn) EQ 0 then epsfcn = MACHEP
  if n_elements(nfev) EQ 0 then nfev = 0L
  if n_elements(xall) EQ 0 then xall = x
  if n_elements(ifree) EQ 0 then ifree = lindgen(n_elements(xall))
  if n_elements(step) EQ 0 then step = x * 0 

  eps = sqrt(max([epsfcn, MACHEP]));
  n = n_elements(x)
  m = n_elements(fvec)
  fjac = reform(replicate(fvec(0), m, n), m, n, /overwrite)

  h = eps * abs(x)
  wh = where(h EQ 0, ct)
  if ct GT 0 then h(wh) = eps

  ;; if STEP is given, use that
  ct = 0 & if n_elements(step) GT 0 then wh = where(step GT 0, ct)
  if ct GT 0 then h(wh) = step(wh)

  ;; if LIMITS specified, then respect those
  ct = 0 & if n_elements(ulimited) GT 0 AND n_elements(ulimit) GT 0 then $
    wh = where(ulimited AND (x GT ulimit-h), ct)
  if ct GT 0 then h(wh) = -h(wh)

  for j=0L, n-1 do begin
      xp = xall
      xp(ifree(j)) = xp(ifree(j)) + h(j)
      
      !err = 0
      fp = call_function(fcn, xp, _EXTRA=fcnargs)
      
      nfev = nfev + 1
      iflag = !err
      !err = 0
      if iflag LT 0 then return, !values.d_nan

      fjac(*,j) = (fp-fvec)/h(j)
  endfor

  return, fjac
end

;     **********
;
;     function enorm
;
;     given an n-vector x, this function calculates the
;     euclidean norm of x.
;
;     the euclidean norm is computed by accumulating the sum of
;     squares in three different sums. the sums of squares for the
;     small and large components are scaled so that no overflows
;     occur. non-destructive underflows are permitted. underflows
;     and overflows do not occur in the computation of the unscaled
;     sum of squares for the intermediate components.
;     the definitions of small, intermediate and large components
;     depend on two ;constants, rdwarf and rgiant. the main
;     restrictions on these constants are that rdwarf**2 not
;     underflow and rgiant**2 not overflow. the constants
;     given here are suitable for every known computer.
;
;     the function statement is
;
;	double precision function enorm(n,x)
;
;     where
;
;	n is a positive integer input variable.
;
;	x is an input array of length n.
;
;     subprograms called
;
;	fortran-supplied ... dabs,dsqrt
;
;     argonne national laboratory. minpack project. march 1980.
;     burton s. garbow, kenneth e. hillstrom, jorge j. more
;
function mpfit_enorm, vec, tricky=tricky
; Very simple-minded sum-of-squares
;  ans0 = sqrt(total(vec^2, 1))

  sz = size(vec)
  mch = machar(double=(sz(sz(0)+1) EQ 5))
  MACHEP = mch.eps
  rdwarf = sqrt(mch.xmin) * 0.9
  rgiant = sqrt(mch.xmax) * 0.9
  sz = size(vec)
  if sz(0) EQ 0 then return, abs(vec)
  if sz(0) EQ 1 then begin
      nr = 1L
      nc = sz(1)
  endif
  if sz(0) EQ 2 then begin
      nr = sz(2)
      nc = sz(1)
  endif
  if sz(0) EQ 2 AND (sz(1) EQ 1) then begin
      vec = vec(*)
      nr = 1L
      nc = n_elements(vec)
  endif
  ans = replicate(vec(0)*0, nr)
  zero = vec(0)*0
  if n_elements(ans) EQ 1 then ans = zero

  for j = 0, nr-1 do begin
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      agiant = rgiant/float(nc)
      x = vec(*,j)
      xabs = abs(x)

      big = (xabs GE agiant)
      sml = (xabs LE rdwarf AND xabs GT 0)

      wh = where(NOT big AND NOT sml, ct)
      if ct GT 0 then s2 = total(xabs(wh)^2)

      wh = where(big, ct)
      if ct GT 0 then begin
          x1max = max(xabs(wh))
          s1 = total((xabs(wh)/x1max)^2)
      endif

      wh = where(sml, ct)
      if ct GT 0 then begin
          x3max = max(xabs(wh))
          s3 = total((xabs(wh)/x3max)^2)
      endif

      if s1 NE 0 then begin
          ans(j) = x1max*sqrt(s1 + (s2/x1max)/x1max)
      endif else if s2 NE 0 then begin
          if s2 GE x3max then $
            temp = s2*((x3max/s2)*(x3max*s3)+1) $
          else $
            temp = x3max*((s2/x3max)+(x3max*s3))
          ans(j) = sqrt(temp)
      endif else begin
          ans(j) = x3max*sqrt(s3)
      endelse
  endfor

  return, ans
end

;     **********
;
;     subroutine qrfac
;
;     this subroutine uses householder transformations with column
;     pivoting (optional) to compute a qr factorization of the
;     m by n matrix a. that is, qrfac determines an orthogonal
;     matrix q, a permutation matrix p, and an upper trapezoidal
;     matrix r with diagonal elements of nonincreasing magnitude,
;     such that a*p = q*r. the householder transformation for
;     column k, k = 1,2,...,min(m,n), is of the form
;
;			    t
;	    i - (1/u(k))*u*u
;
;     where u has zeros in the first k-1 positions. the form of
;     this transformation and the method of pivoting first
;     appeared in the corresponding linpack subroutine.
;
;     the subroutine statement is
;
;	subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
;
;     where
;
;	m is a positive integer input variable set to the number
;	  of rows of a.
;
;	n is a positive integer input variable set to the number
;	  of columns of a.
;
;	a is an m by n array. on input a contains the matrix for
;	  which the qr factorization is to be computed. on output
;	  the strict upper trapezoidal part of a contains the strict
;	  upper trapezoidal part of r, and the lower trapezoidal
;	  part of a contains a factored form of q (the non-trivial
;	  elements of the u vectors described above).
;
;	lda is a positive integer input variable not less than m
;	  which specifies the leading dimension of the array a.
;
;	pivot is a logical input variable. if pivot is set true,
;	  then column pivoting is enforced. if pivot is set false,
;	  then no column pivoting is done.
;
;	ipvt is an integer output array of length lipvt. ipvt
;	  defines the permutation matrix p such that a*p = q*r.
;	  column j of p is column ipvt(j) of the identity matrix.
;	  if pivot is false, ipvt is not referenced.
;
;	lipvt is a positive integer input variable. if pivot is false,
;	  then lipvt may be as small as 1. if pivot is true, then
;	  lipvt must be at least n.
;
;	rdiag is an output array of length n which contains the
;	  diagonal elements of r.
;
;	acnorm is an output array of length n which contains the
;	  norms of the corresponding columns of the input matrix a.
;	  if this information is not needed, then acnorm can coincide
;	  with rdiag.
;
;	wa is a work array of length n. if pivot is false, then wa
;	  can coincide with rdiag.
;
;     subprograms called
;
;	minpack-supplied ... dpmpar,enorm
;
;	fortran-supplied ... dmax1,dsqrt,min0
;
;     argonne national laboratory. minpack project. march 1980.
;     burton s. garbow, kenneth e. hillstrom, jorge j. more
;
;     **********
pro mpfit_qrfac, a, ipvt, rdiag, acnorm, pivot=pivot

  sz = size(a)
  mch = machar(double=(sz(sz(0)+1) EQ 5))
  MACHEP = mch.eps
  DWARF = mch.xmin
  sz = size(a)
  m = sz(1)
  n = sz(2)
  
  ;; Compute the initial column norms and initialize arrays
  acnorm = mpfit_enorm(a)
  rdiag = acnorm
  wa = rdiag
  if keyword_set(pivot) then ipvt = lindgen(n)
  
  ;; Reduce a to r with householder transformations
  minmn = min([m,n])
  for j = 0L, minmn-1 do begin
      if NOT keyword_set(pivot) then goto, HOUSE1
      
      ;; Bring the column of largest norm into the pivot position
      rmax = max(rdiag(j:*))
      kmax = where(rdiag(j:*) EQ rmax, ct) + j
      if ct LE 0 then goto, HOUSE1
      kmax = kmax(0)
      
      ;; Exchange rows
      temp       = a(*,j)    & a(*,j)     = a(*,kmax)  & a(*,kmax)  = temp
      rdiag(kmax) = rdiag(j)
      wa(kmax)   = wa(j)
      temp       = ipvt(j)   & ipvt(j)    = ipvt(kmax) & ipvt(kmax) = temp
      
      HOUSE1:

      ;; Compute the householder transformation to reduce the jth
      ;; column of a to a multiple of the jth unit vector
      ajnorm = mpfit_enorm(a(j:*,j))
      if ajnorm EQ 0 then goto, NEXT_ROW
      if a(j,j) LT 0 then ajnorm = -ajnorm
      
      a(j:*,j) = a(j:*,j) / ajnorm
      a(j,j) = a(j,j) + 1
      
      ;; Apply the transformation to the remaining columns
      ;; and update the norms
      if j+1 LT n then begin
          for k=j+1, n-1 do begin
              sum = total(a(j:*,k)*a(j:*,j))
              temp = sum/a(j,j)
              a(j:*,k) = a(j:*,k) - temp * a(j:*,j)

              if keyword_set(pivot) AND rdiag(k) NE 0 then begin
                  temp = a(j,k)/rdiag(k)
                  rdiag(k) = rdiag(k) * sqrt(max([-temp*temp + 0, 1]))
                  temp = rdiag(k)/wa(k)
                  if 0.05D*temp*temp LE MACHEP then begin
                      rdiag(k) = mpfit_enorm(a((j+1):*,k))
                      wa(k) = rdiag(k)
                  endif
              endif
          endfor
      endif

      NEXT_ROW:
      rdiag(j) = -ajnorm
  endfor

  return
end

;     **********
;
;     subroutine qrsolv
;
;     given an m by n matrix a, an n by n diagonal matrix d,
;     and an m-vector b, the problem is to determine an x which
;     solves the system
;
;           a*x = b ,     d*x = 0 ,
;
;     in the least squares sense.
;
;     this subroutine completes the solution of the problem
;     if it is provided with the necessary information from the
;     qr factorization, with column pivoting, of a. that is, if
;     a*p = q*r, where p is a permutation matrix, q has orthogonal
;     columns, and r is an upper triangular matrix with diagonal
;     elements of nonincreasing magnitude, then qrsolv expects
;     the full upper triangle of r, the permutation matrix p,
;     and the first n components of (q transpose)*b. the system
;     a*x = b, d*x = 0, is then equivalent to
;
;                  t       t
;           r*z = q *b ,  p *d*p*z = 0 ,
;
;     where x = p*z. if this system does not have full rank,
;     then a least squares solution is obtained. on output qrsolv
;     also provides an upper triangular matrix s such that
;
;            t   t               t
;           p *(a *a + d*d)*p = s *s .
;
;     s is computed within qrsolv and may be of separate interest.
;
;     the subroutine statement is
;
;       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
;
;     where
;
;       n is a positive integer input variable set to the order of r.
;
;       r is an n by n array. on input the full upper triangle
;         must contain the full upper triangle of the matrix r.
;         on output the full upper triangle is unaltered, and the
;         strict lower triangle contains the strict upper triangle
;         (transposed) of the upper triangular matrix s.
;
;       ldr is a positive integer input variable not less than n
;         which specifies the leading dimension of the array r.
;
;       ipvt is an integer input array of length n which defines the
;         permutation matrix p such that a*p = q*r. column j of p
;         is column ipvt(j) of the identity matrix.
;
;       diag is an input array of length n which must contain the
;         diagonal elements of the matrix d.
;
;       qtb is an input array of length n which must contain the first
;         n elements of the vector (q transpose)*b.
;
;       x is an output array of length n which contains the least
;         squares solution of the system a*x = b, d*x = 0.
;
;       sdiag is an output array of length n which contains the
;         diagonal elements of the upper triangular matrix s.
;
;       wa is a work array of length n.
;
;     subprograms called
;
;       fortran-supplied ... dabs,dsqrt
;
;     argonne national laboratory. minpack project. march 1980.
;     burton s. garbow, kenneth e. hillstrom, jorge j. more
;
pro mpfit_qrsolv, r, ipvt, diag, qtb, x, sdiag

  sz = size(r)
  m = sz(1)
  n = sz(2)
  delm = lindgen(n) * (m+1) ;; Diagonal elements of r

  ;; copy r and (q transpose)*b to preserve input and initialize s.
  ;; in particular, save the diagonal elements of r in x.

  for j = 0L, n-1 do $
    r(j:n-1,j) = r(j,j:n-1)
  x = r(delm)
  wa = qtb
  zero = qtb(0)*0
  half = zero + 0.5
  quart = zero + 0.25

  ;; Eliminate the diagonal matrix d using a givens rotation
  for j = 0L, n-1 do begin
      l = ipvt(j)
      if diag(l) EQ 0 then goto, STORE_RESTORE
      sdiag(j:*) = 0
      sdiag(j) = diag(l)

      ;; The transformations to eliminate the row of d modify only a
      ;; single element of (q transpose)*b beyond the first n, which
      ;; is initially zero.

      qtbpj = zero
      for k = j, n-1 do begin
          if sdiag(k) EQ 0 then goto, ELIM_NEXT_LOOP
          if abs(r(k,k)) LT abs(sdiag(k)) then begin
              cotan = r(k,k)/sdiag(k)
              sin = half/sqrt(quart + quart*cotan*cotan)
              cos = sin*cotan
          endif else begin
              tan = sdiag(k)/r(k,k)
              cos = half/sqrt(quart + quart*tan*tan)
              sin = cos*tan
          endelse
          
          ;; Compute the modified diagonal element of r and the
          ;; modified element of ((q transpose)*b,0).
          r(k,k) = cos*r(k,k) + sin*sdiag(k)
          temp = cos*wa(k) + sin*qtbpj
          qtbpj = -sin*wa(k) + cos*qtbpj
          wa(k) = temp

          ;; Accumulate the transformation in the row of s
          if n GT k+1 then begin
              temp = cos*r(k+1:n-1,k) + sin*sdiag(k+1:n-1)
              sdiag(k+1:n-1) = -sin*r(k+1:n-1,k) + cos*sdiag(k+1:n-1)
              r(k+1:n-1,k) = temp
          endif
ELIM_NEXT_LOOP:
      endfor

STORE_RESTORE:
      sdiag(j) = r(j,j)
      r(j,j) = x(j)
  endfor

  ;; Solve the triangular system for z.  If the system is singular
  ;; then obtain a least squares solution
  nsing = n
  wh = where(sdiag EQ 0, ct)
  if ct GT 0 then begin
      nsing = wh(0)
      wa(nsing:*) = 0
  endif

  if nsing GE 1 then begin
      wa(nsing-1) = wa(nsing-1)/sdiag(nsing-1) ;; Degenerate case
      for j=nsing-2,0,-1 do begin  ;; Reverse loop
          sum = total(r(j+1:nsing-1,j)*wa(j+1:nsing-1))
          wa(j) = (wa(j)-sum)/sdiag(j)
      endfor
  endif

  ;; Permute the components of z back to components of x
  x(ipvt) = wa

  return
end
      
  
;
;     subroutine lmpar
;
;     given an m by n matrix a, an n by n nonsingular diagonal
;     matrix d, an m-vector b, and a positive number delta,
;     the problem is to determine a value for the parameter
;     par such that if x solves the system
;
;	    a*x = b ,	  sqrt(par)*d*x = 0 ,
;
;     in the least squares sense, and dxnorm is the euclidean
;     norm of d*x, then either par is zero and
;
;	    (dxnorm-delta) .le. 0.1*delta ,
;
;     or par is positive and
;
;	    abs(dxnorm-delta) .le. 0.1*delta .
;
;     this subroutine completes the solution of the problem
;     if it is provided with the necessary information from the
;     qr factorization, with column pivoting, of a. that is, if
;     a*p = q*r, where p is a permutation matrix, q has orthogonal
;     columns, and r is an upper triangular matrix with diagonal
;     elements of nonincreasing magnitude, then lmpar expects
;     the full upper triangle of r, the permutation matrix p,
;     and the first n components of (q transpose)*b. on output
;     lmpar also provides an upper triangular matrix s such that
;
;	     t	 t		     t
;	    p *(a *a + par*d*d)*p = s *s .
;
;     s is employed within lmpar and may be of separate interest.
;
;     only a few iterations are generally needed for convergence
;     of the algorithm. if, however, the limit of 10 iterations
;     is reached, then the output par will contain the best
;     value obtained so far.
;
;     the subroutine statement is
;
;	subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
;			 wa1,wa2)
;
;     where
;
;	n is a positive integer input variable set to the order of r.
;
;	r is an n by n array. on input the full upper triangle
;	  must contain the full upper triangle of the matrix r.
;	  on output the full upper triangle is unaltered, and the
;	  strict lower triangle contains the strict upper triangle
;	  (transposed) of the upper triangular matrix s.
;
;	ldr is a positive integer input variable not less than n
;	  which specifies the leading dimension of the array r.
;
;	ipvt is an integer input array of length n which defines the
;	  permutation matrix p such that a*p = q*r. column j of p
;	  is column ipvt(j) of the identity matrix.
;
;	diag is an input array of length n which must contain the
;	  diagonal elements of the matrix d.
;
;	qtb is an input array of length n which must contain the first
;	  n elements of the vector (q transpose)*b.
;
;	delta is a positive input variable which specifies an upper
;	  bound on the euclidean norm of d*x.
;
;	par is a nonnegative variable. on input par contains an
;	  initial estimate of the levenberg-marquardt parameter.
;	  on output par contains the final estimate.
;
;	x is an output array of length n which contains the least
;	  squares solution of the system a*x = b, sqrt(par)*d*x = 0,
;	  for the output par.
;
;	sdiag is an output array of length n which contains the
;	  diagonal elements of the upper triangular matrix s.
;
;	wa1 and wa2 are work arrays of length n.
;
;     subprograms called
;
;	minpack-supplied ... dpmpar,enorm,qrsolv
;
;	fortran-supplied ... dabs,dmax1,dmin1,dsqrt
;
;     argonne national laboratory. minpack project. march 1980.
;     burton s. garbow, kenneth e. hillstrom, jorge j. more
;
function mpfit_lmpar, r, ipvt, diag, qtb, delta, x, sdiag, par=par

  sz = size(r)
  mch = machar(double=(sz(sz(0)+1) EQ 5))
  MACHEP = mch.eps
  DWARF = mch.xmin
  sz = size(r)
  m = sz(1)
  n = sz(2)
  delm = lindgen(n) * (m+1) ;; Diagonal elements of r

  ;; Compute and store in x the gauss-newton direction.  If the
  ;; jacobian is rank-deficient, obtain a least-squares solution
  nsing = n
  wa1 = qtb
  wh = where(r(delm) EQ 0, ct)
  if ct GT 0 then begin
      nsing = wh(0)
      wa1(nsing:*) = 0
  endif

  if nsing GE 1 then begin
      for j=nsing-1,0,-1 do begin  ;; Reverse loop
          wa1(j) = wa1(j)/r(j,j)
          if (j-1 GE 0) then $
            wa1(0:(j-1)) = wa1(0:(j-1)) - r(0:(j-1),j)*wa1(j)
      endfor
  endif

  ;; Note: ipvt here is a permutation array
  x(ipvt) = wa1

  ;; Initialize the iteration counter.  Evaluate the function at the
  ;; origin, and test for acceptance of the gauss-newton direction
  iter = 0L
  wa2 = diag * x
  dxnorm = mpfit_enorm(wa2)
  fp = dxnorm - delta
  if fp LE 0.1 then goto, TERMINATE

  ;; If the jacobian is not rank deficient, the newton step provides a
  ;; lower bound, parl, for the zero of the function.  Otherwise set
  ;; this bound to zero.
  
  zero = wa2(0)*0
  parl = zero
  if nsing GE n then begin
      wa1 = diag(ipvt)*wa2(ipvt)/dxnorm

      wa1(0) = wa1(0) / r(0,0) ;; Degenerate case 
      for j=1L, n-1 do begin   ;; Note "1" here, not zero
          sum = total(r(0:(j-1),j)*wa1(0:(j-1)))
          wa1(j) = (wa1(j) - sum)/r(j,j)
      endfor

      temp = mpfit_enorm(wa1)
      parl = ((fp/delta)/temp)/temp
  endif

  ;; Calculate an upper bound, paru, for the zero of the function
  for j=0, n-1 do begin
      sum = total(r(0:j,j)*qtb(0:j))
      wa1(j) = sum/diag(ipvt(j))
  endfor
  gnorm = mpfit_enorm(wa1)
  paru  = gnorm/delta
  if paru EQ 0 then paru = DWARF/min([delta,0.1])

  ;; If the input par lies outside of the interval (parl,paru), set
  ;; par to the closer endpoint

  par = max([par,parl])
  par = min([par,paru])
  if par EQ 0 then par = gnorm/dxnorm

  ;; Beginning of an interation
  ITERATION:
  iter = iter + 1
  
  ;; Evaluate the function at the current value of par
  if par EQ 0 then par = max([DWARF, paru*0.001])
  temp = sqrt(par)
  wa1 = temp * diag
  mpfit_qrsolv, r, ipvt, wa1, qtb, x, sdiag
  wa2 = diag*x
  dxnorm = mpfit_enorm(wa2)
  temp = fp
  fp = dxnorm - delta

  if (abs(fp) LE 0.1D*delta) $
    OR ((parl EQ 0) AND (fp LE temp) AND (temp LT 0)) $
    OR (iter EQ 10) then goto, TERMINATE

  ;; Compute the newton correction
  wa1 = diag(ipvt)*wa2(ipvt)/dxnorm

  for j=0,n-2 do begin
      wa1(j) = wa1(j)/sdiag(j)
      wa1(j+1:n-1) = wa1(j+1:n-1) - r(j+1:n-1,j)*wa1(j)
  endfor
  wa1(n-1) = wa1(n-1)/sdiag(n-1) ;; Degenerate case

  temp = mpfit_enorm(wa1)
  parc = ((fp/delta)/temp)/temp

  ;; Depending on the sign of the function, update parl or paru
  if fp GT 0 then parl = max([parl,par])
  if fp LT 0 then paru = min([paru,par])

  ;; Compute an improved estimate for par
  par = max([parl, par+parc])

  ;; End of an iteration
  goto, ITERATION
  
TERMINATE:
  ;; Termination
  if iter EQ 0 then return, par(0)*0
  return, par
end

pro mpfit_defiter, fcn, x, iter, FUNCTARGS=fcnargs, fmt=fmt, $
         quiet=quiet, _EXTRA=iterargs

  if keyword_set(quiet) then return
  fvec = call_function(fcn, x, _EXTRA=fcnargs)
  fnorm = mpfit_enorm(fvec)

  print, iter, fnorm^2, $
    format='("Iter ",I6,"   CHI-SQUARE = ",G20.8)'
  if n_elements(fmt) GT 0 then begin
      print, p, format=fmt
  endif else begin
      p = '  P('+strtrim(lindgen(n_elements(x)),2)+') = ' $
        + strtrim(string(x,format='(G20.6)'),2) + '  '
      print, '         '+p, format='(A)'
  endelse
  
  return
end

;     **********
;
;     subroutine covar
;
;     given an m by n matrix a, the problem is to determine
;     the covariance matrix corresponding to a, defined as
;
;                    t
;           inverse(a *a) .
;
;     this subroutine completes the solution of the problem
;     if it is provided with the necessary information from the
;     qr factorization, with column pivoting, of a. that is, if
;     a*p = q*r, where p is a permutation matrix, q has orthogonal
;     columns, and r is an upper triangular matrix with diagonal
;     elements of nonincreasing magnitude, then covar expects
;     the full upper triangle of r and the permutation matrix p.
;     the covariance matrix is then computed as
;
;                      t     t
;           p*inverse(r *r)*p  .
;
;     if a is nearly rank deficient, it may be desirable to compute
;     the covariance matrix corresponding to the linearly independent
;     columns of a. to define the numerical rank of a, covar uses
;     the tolerance tol. if l is the largest integer such that
;
;           abs(r(l,l)) .gt. tol*abs(r(1,1)) ,
;
;     then covar computes the covariance matrix corresponding to
;     the first l columns of r. for k greater than l, column
;     and row ipvt(k) of the covariance matrix are set to zero.
;
;     the subroutine statement is
;
;       subroutine covar(n,r,ldr,ipvt,tol,wa)
;
;     where
;
;       n is a positive integer input variable set to the order of r.
;
;       r is an n by n array. on input the full upper triangle must
;         contain the full upper triangle of the matrix r. on output
;         r contains the square symmetric covariance matrix.
;
;       ldr is a positive integer input variable not less than n
;         which specifies the leading dimension of the array r.
;
;       ipvt is an integer input array of length n which defines the
;         permutation matrix p such that a*p = q*r. column j of p
;         is column ipvt(j) of the identity matrix.
;
;       tol is a nonnegative input variable used to define the
;         numerical rank of a in the manner described above.
;
;       wa is a work array of length n.
;
;     subprograms called
;
;       fortran-supplied ... dabs
;
;     argonne national laboratory. minpack project. august 1980.
;     burton s. garbow, kenneth e. hillstrom, jorge j. more
;
;     **********
function mpfit_covar, rr, ipvt, tol=tol

  if n_elements(tol) EQ 0 then tol = 1.D-14
  sz = size(rr)
  if sz(0) NE 2 then begin
      message, 'ERROR: r must be a two-dimensional matrix'
      return, -1L
  endif
  n = sz(1)
  if n NE sz(2) then begin
      message, 'ERROR: r must be a square matrix'
      return, -1L
  endif

  if n_elements(ipvt) EQ 0 then ipvt = lindgen(n)
  r = rr
  
  ;; For the inverse of r in the full upper triangle of r
  l = -1L
  tolr = tol * abs(r(0,0))
  zero = r(0,0) * 0.
  one  = zero + 1.
  for k = 0L, n-1 do begin
      if abs(r(k,k)) LE tolr then goto, INV_END_LOOP
      r(k,k) = one/r(k,k)
      for j = 0L, k-1 do begin
          temp = r(k,k) * r(j,k)
          r(j,k) = zero
          r(0:j,k) = r(0:j,k) - temp*r(0:j,j)
      endfor
      l = k
  endfor
  INV_END_LOOP:

  ;; Form the full upper triangle of the inverse of (r transpose)*r
  ;; in the full upper triangle of r
  if l GE 0 then $
    for k = 0L, l do begin
      for j = 0L, k-1 do begin
          temp = r(j,k)
          r(0:j,j) = r(0:j,j) + temp*r(0:j,k)
      endfor
      temp = r(k,k)
      r(0:k,k) = temp * r(0:k,k)
  endfor

  ;; For the full lower triangle of the covariance matrix
  ;; in the strict lower triangle or and in wa
  wa = replicate(r(0,0), n)
  for j = 0L, n-1 do begin
      jj = ipvt(j)
      sing = j GT l
      for i = 0L, j do begin
          if sing then r(i,j) = zero
          ii = ipvt(i)
          if ii GT jj then r(ii,jj) = r(i,j)
          if ii LT jj then r(jj,ii) = r(i,j)
      endfor
      wa(jj) = r(j,j)
  endfor

  ;; Symmetrize the covariance matrix in r
  for j = 0L, n-1 do begin
      r(0:j,j) = r(j,0:j)
      r(j,j) = wa(j)
  endfor

  return, r
end

;     **********
;
;     subroutine lmdif
;
;     the purpose of lmdif is to minimize the sum of the squares of
;     m nonlinear functions in n variables by a modification of
;     the levenberg-marquardt algorithm. the user must provide a
;     subroutine which calculates the functions. the jacobian is
;     then calculated by a forward-difference approximation.
;
;     the subroutine statement is
;
;	subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
;			 diag,mode,factor,nprint,info,nfev,fjac,
;			 ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
;
;     where
;
;	fcn is the name of the user-supplied subroutine which
;	  calculates the functions. fcn must be declared
;	  in an external statement in the user calling
;	  program, and should be written as follows.
;
;	  subroutine fcn(m,n,x,fvec,iflag)
;	  integer m,n,iflag
;	  double precision x(n),fvec(m)
;	  ----------
;	  calculate the functions at x and
;	  return this vector in fvec.
;	  ----------
;	  return
;	  end
;
;	  the value of iflag should not be changed by fcn unless
;	  the user wants to terminate execution of lmdif.
;	  in this case set iflag to a negative integer.
;
;	m is a positive integer input variable set to the number
;	  of functions.
;
;	n is a positive integer input variable set to the number
;	  of variables. n must not exceed m.
;
;	x is an array of length n. on input x must contain
;	  an initial estimate of the solution vector. on output x
;	  contains the final estimate of the solution vector.
;
;	fvec is an output array of length m which contains
;	  the functions evaluated at the output x.
;
;	ftol is a nonnegative input variable. termination
;	  occurs when both the actual and predicted relative
;	  reductions in the sum of squares are at most ftol.
;	  therefore, ftol measures the relative error desired
;	  in the sum of squares.
;
;	xtol is a nonnegative input variable. termination
;	  occurs when the relative error between two consecutive
;	  iterates is at most xtol. therefore, xtol measures the
;	  relative error desired in the approximate solution.
;
;	gtol is a nonnegative input variable. termination
;	  occurs when the cosine of the angle between fvec and
;	  any column of the jacobian is at most gtol in absolute
;	  value. therefore, gtol measures the orthogonality
;	  desired between the function vector and the columns
;	  of the jacobian.
;
;	maxfev is a positive integer input variable. termination
;	  occurs when the number of calls to fcn is at least
;	  maxfev by the end of an iteration.
;
;	epsfcn is an input variable used in determining a suitable
;	  step length for the forward-difference approximation. this
;	  approximation assumes that the relative errors in the
;	  functions are of the order of epsfcn. if epsfcn is less
;	  than the machine precision, it is assumed that the relative
;	  errors in the functions are of the order of the machine
;	  precision.
;
;	diag is an array of length n. if mode = 1 (see
;	  below), diag is internally set. if mode = 2, diag
;	  must contain positive entries that serve as
;	  multiplicative scale factors for the variables.
;
;	mode is an integer input variable. if mode = 1, the
;	  variables will be scaled internally. if mode = 2,
;	  the scaling is specified by the input diag. other
;	  values of mode are equivalent to mode = 1.
;
;	factor is a positive input variable used in determining the
;	  initial step bound. this bound is set to the product of
;	  factor and the euclidean norm of diag*x if nonzero, or else
;	  to factor itself. in most cases factor should lie in the
;	  interval (.1,100.). 100. is a generally recommended value.
;
;	nprint is an integer input variable that enables controlled
;	  printing of iterates if it is positive. in this case,
;	  fcn is called with iflag = 0 at the beginning of the first
;	  iteration and every nprint iterations thereafter and
;	  immediately prior to return, with x and fvec available
;	  for printing. if nprint is not positive, no special calls
;	  of fcn with iflag = 0 are made.
;
;	info is an integer output variable. if the user has
;	  terminated execution, info is set to the (negative)
;	  value of iflag. see description of fcn. otherwise,
;	  info is set as follows.
;
;	  info = 0  improper input parameters.
;
;	  info = 1  both actual and predicted relative reductions
;		    in the sum of squares are at most ftol.
;
;	  info = 2  relative error between two consecutive iterates
;		    is at most xtol.
;
;	  info = 3  conditions for info = 1 and info = 2 both hold.
;
;	  info = 4  the cosine of the angle between fvec and any
;		    column of the jacobian is at most gtol in
;		    absolute value.
;
;	  info = 5  number of calls to fcn has reached or
;		    exceeded maxfev.
;
;	  info = 6  ftol is too small. no further reduction in
;		    the sum of squares is possible.
;
;	  info = 7  xtol is too small. no further improvement in
;		    the approximate solution x is possible.
;
;	  info = 8  gtol is too small. fvec is orthogonal to the
;		    columns of the jacobian to machine precision.
;
;	nfev is an integer output variable set to the number of
;	  calls to fcn.
;
;	fjac is an output m by n array. the upper n by n submatrix
;	  of fjac contains an upper triangular matrix r with
;	  diagonal elements of nonincreasing magnitude such that
;
;		 t     t	   t
;		p *(jac *jac)*p = r *r,
;
;	  where p is a permutation matrix and jac is the final
;	  calculated jacobian. column j of p is column ipvt(j)
;	  (see below) of the identity matrix. the lower trapezoidal
;	  part of fjac contains information generated during
;	  the computation of r.
;
;	ldfjac is a positive integer input variable not less than m
;	  which specifies the leading dimension of the array fjac.
;
;	ipvt is an integer output array of length n. ipvt
;	  defines a permutation matrix p such that jac*p = q*r,
;	  where jac is the final calculated jacobian, q is
;	  orthogonal (not stored), and r is upper triangular
;	  with diagonal elements of nonincreasing magnitude.
;	  column j of p is column ipvt(j) of the identity matrix.
;
;	qtf is an output array of length n which contains
;	  the first n elements of the vector (q transpose)*fvec.
;
;	wa1, wa2, and wa3 are work arrays of length n.
;
;	wa4 is a work array of length m.
;
;     subprograms called
;
;	user-supplied ...... fcn
;
;	minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
;
;	fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
;
;     argonne national laboratory. minpack project. march 1980.
;     burton s. garbow, kenneth e. hillstrom, jorge j. more
;
;     **********
function mpfit, fcn, xall, FUNCTARGS=fcnargs, $
                ftol=ftol, xtol=xtol, gtol=gtol, epsfcn=epsfcn, $
                nfev=nfev, maxiter=maxiter, errmsg=errmsg, $
                factor=factor, nprint=nprint, STATUS=info, $
                iterproc=iterproc, iterargs=iterargs, $
                diag=diag, rescale=rescale, $
                perror=perror, covar=covar, nocovar=nocovar, bestnorm=fnorm, $
                parinfo=parinfo, quiet=quiet

  if n_elements(ftol) EQ 0 then ftol = 1.D-10
  if n_elements(xtol) EQ 0 then xtol = 1.D-10
  if n_elements(gtol) EQ 0 then gtol = 1.D-10
  if n_elements(factor) EQ 0 then factor = 100.
  if n_elements(nprint) EQ 0 then nprint = 1
  if n_elements(iterproc) EQ 0 then iterproc = 'mpfit_defiter'

  info = 0L
  iflag = 0L
  nfev = 0L
  errmsg = ''
  fnorm  = -1.D
  fnorm1 = -1.D

  ;; Handle error conditions gracefully
  catch, catcherror
  if catcherror NE 0 then begin
      catch, /cancel
      message, !err_string, /info
      message, 'Error condition detected. Returning to MAIN level.', /info
      return, !values.d_nan
  endif

  ;; Parinfo:
  ;; --------------- STARTING/CONFIG INFO (passed in to routine, not changed)
  ;; .value   - starting value for parameter
  ;; .fixed   - parameter is fixed
  ;; .limited - a two-element array, if parameter is bounded on
  ;;            lower/upper side
  ;; .limits - a two-element array, lower/upper parameter bounds, if
  ;;           limited vale is set
  ;; .step   - step size in Jacobian calc, if greater than zero

  ;; Parameters can either be stored in parinfo, or x.  Parinfo takes
  ;; precedence if it exists.
  if n_elements(xall) EQ 0 AND n_elements(parinfo) EQ 0 then begin
      errmsg = 'ERROR: must pass parameters in P or PARINFO'
      goto, TERMINATE
  endif
  if n_elements(xall) EQ 0 then begin
      xall = parinfo(*).value
      sz = size(xall)
      ;; Convert to double if parameters are not float or double
      if sz(sz(0)+1) NE 4 AND sz(sz(0)+1) NE 5 then $
        xall = double(xall)
  endif
  if n_elements(parinfo) EQ 0 then begin
      parinfo = replicate({value:0.D, fixed:0, $
                           limited:[0,0], limits:[0.D,0], step:0.D}, $
                           n_elements(xall))
      parinfo(*).value = xall
  endif
  ;; Error checking on parinfo
  wh = where( (parinfo(*).limited(0) AND xall LT parinfo(*).limits(0)) $
              OR (parinfo(*).limited(1) AND xall GT parinfo(*).limits(1)), ct)
  if ct GT 0 then begin
      errmsg = 'ERROR: parameters are not within PARINFO limits'
      goto, TERMINATE
  endif
  wh = where(parinfo(*).limited(0) AND parinfo(*).limited(1) $
             AND parinfo(*).limits(0) GE parinfo(*).limits(1) $
             AND NOT parinfo(*).fixed, ct)
  if ct GT 0 then begin
      errmsg = 'ERROR: PARINFO parameter limits are not consistent'
      goto, TERMINATE
  endif

  ;; Get freely adjustable parameters
  ifree = where(parinfo(*).fixed EQ 0, ct)
  if ct EQ 0 then begin
      errmsg = 'ERROR: no free parameters'
      goto, TERMINATE
  endif
  xnew = xall      ;; xnew is the set of parameters to be returned
  x = xnew(ifree)  ;; x is the set of free parameters
  qulim = parinfo(ifree).limited(1)
  ulim  = parinfo(ifree).limits(1)
  qllim = parinfo(ifree).limited(0)
  llim  = parinfo(ifree).limits(0)
  step  = parinfo(ifree).step

  n = n_elements(x)
  if n_elements(maxiter) EQ 0 then maxiter = 200L

  ;; Check input parameters for errors
  if (n LE 0) OR (ftol LT 0) OR (xtol LT 0) OR (gtol LT 0) $
    OR (maxiter LE 0) OR (factor LE 0) then begin
      errmsg = 'ERROR: input keywords are inconsistent'
      goto, TERMINATE
  endif
 
  if keyword_set(rescale) then begin
      errmsg = 'ERROR: DIAG parameter scales are inconsistent'
      if n_elements(diag) LT n then goto, TERMINATE
      wh = where(diag LE 0, ct)
      if ct GT 0 then goto, TERMINATE
      errmsg = ''
  endif

  !err = 1
  fvec = call_function(fcn, xnew, _EXTRA=fcnargs)
  nfev = 1
  iflag = !err
  !err = 0
  if iflag LT 0 then begin
      errmsg = 'ERROR: first call to "'+fcn+'" failed'
      goto, TERMINATE
  endif

  sz = size(fvec(0))
  mch = machar(double=(sz(sz(0)+1) EQ 5))
  MACHEP = mch.eps
  DWARF = mch.xmin
  szx = size(x)
  ;; The parameters and the squared deviations should have the same
  ;; type.  Otherwise the MACHAR-based evaluation will fail.
  if sz(sz(0)+1) EQ 5 AND szx(szx(0)+1) NE 5 then begin
      if NOT keyword_set(quiet) then begin
          message, 'WARNING: data is DOUBLE but parameters are FLOAT', /info
          message, '         (converting parameters to DOUBLE)', /info
      endif
      x = double(x)
      xnew = double(xnew)
  endif

  m = n_elements(fvec)
  if (m LT n) then begin
      errmsg = 'ERROR: number of parameters must not exceed data'
      goto, TERMINATE
  endif
  
  fnorm = mpfit_enorm(fvec)

  ;; Initialize Levelberg-Marquardt parameter and iteration counter

  par = x(0) * 0
  iter = 1L
  qtf = x * 0.

  ;; Beginning of the outer loop
  
  OUTER_LOOP:

  ;; If requested, call fcn to enable printing of iterates

  xnew(ifree) = x
  if nprint GT 0 AND iterproc NE '' then begin
      iflag = 0L
      if (iter-1) MOD nprint EQ 0 then begin
          call_procedure, iterproc, fcn, xnew, iter, $
            FUNCTARGS=fcnargs, parinfo=parinfo, quiet=quiet, _EXTRA=iterargs
          iflag = !err
          !err = 0
          if iflag LT 0 then begin
              errmsg = 'WARNING: premature termination by "'+iterproc+'"'
              goto, TERMINATE
          endif
      endif
  endif
  x = xnew(ifree)

  ;; Calculate the jacobian matrix
  iflag = 2
  fjac = mpfit_fdjac2(fcn, x, fvec, step, qulim, ulim, $
                iflag=iflag, epsfcn=epsfcn, nfev=nfev, $
                FUNCTARGS=fcnargs, ifree=ifree, xall=xnew)
  if iflag LT 0 then begin
      errmsg = 'WARNING: premature termination by FDJAC2'
      goto, TERMINATE
  endif

  ;; Determine if any of the parameters are pegged at the limits
  whlpeg = where(qllim AND (x EQ llim), nlpeg)
  whupeg = where(qulim AND (x EQ ulim), nupeg)
  
  ;; See if any "pegged" values should keep their derivatives
  if (nlpeg GT 0) then begin
      ;; Total derivative of sum wrt lower pegged parameters
      for i = 0, nlpeg-1 do begin
          sum = total(fvec * fjac(*,whlpeg(i)))
          if sum GT 0 then fjac(*,whlpeg(i)) = 0
      endfor
  endif
  if (nupeg GT 0) then begin
      ;; Total derivative of sum wrt upper pegged parameters
      for i = 0, nupeg-1 do begin
          sum = total(fvec * fjac(*,whupeg(i)))
          if sum LT 0 then fjac(*,whupeg(i)) = 0
      endfor
  endif

  ;; Compute the QR factorization of the jacobian
  mpfit_qrfac, fjac, ipvt, wa1, wa2, /pivot

  ;; On the first iteration if "diag" is unspecified, scale
  ;; according to the norms of the columns of the initial jacobian
  if (iter EQ 1) then begin

      if NOT keyword_set(rescale) OR (n_elements(diag) LT n) then begin
          diag = wa2
          wh = where (diag EQ 0, ct)
          if ct GT 0 then diag(wh) = 1.D
      endif
      
      ;; On the first iteration, calculate the norm of the scaled x
      ;; and initialize the step bound delta 
      wa3 = diag * x
      xnorm = mpfit_enorm(wa3)
      delta = factor*xnorm
      if delta EQ 0.D then delta = delta(0)*0 + factor
  endif

  ;; Form (q transpose)*fvec and store the first n components in qtf
  wa4 = fvec
  for j=0L, n-1 do begin
      temp3 = fjac(j,j)
      if temp3 NE 0 then begin
          sum = total(fjac(j:*,j)*wa4(j:*))
          temp = -sum/temp3
          wa4(j:*) = wa4(j:*) + fjac(j:*,j) * temp
      endif
      fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)
  endfor

  ;; Compute the norm of the scaled gradient
  gnorm = wa2(0) * 0
  if fnorm NE 0 then begin
      for j=0L, n-1 do begin
          l = ipvt(j)
          if wa2(l) NE 0 then begin
              sum = total(fjac(0:j,j)*qtf(0:j)/fnorm)
              gnorm = max([gnorm,abs(sum/wa2(l))])
          endif
      endfor
  endif

  ;; Test for convergence of the gradient norm
  if gnorm LE gtol then info = 4
  if info NE 0 then goto, TERMINATE

  ;; Rescale if necessary
  if NOT keyword_set(rescale) then $
    diag = diag > wa2

  ;; Beginning of the inner loop
  INNER_LOOP:
  
  ;; Determine the levenberg-marquardt parameter
  par = mpfit_lmpar(fjac, ipvt, diag, qtf, delta, wa1, wa2, par=par)

  ;; Store the direction p and x+p. Calculate the norm of p
  wa1 = -wa1

  ;; Do not allow any steps out of bounds
  if nlpeg GT 0 then wa1(whlpeg) = wa1(whlpeg) > 0
  if nupeg GT 0 then wa1(whupeg) = wa1(whupeg) < 0

  ;; Respect the limits.  If a step were to go out of bounds, then
  ;; we should take a step in the same direction but shorter distance.
  ;; The step should take us right to the limit in that case.
  alpha = x(0)*0 + 1.
  whl = where((abs(wa1) GT MACHEP) AND qllim AND (x + wa1 LT llim), lct)
  if lct GT 0 then $
    alpha = min([alpha, (llim(whl)-x(whl))/wa1(whl)])
  whu = where((abs(wa1) GT MACHEP) AND qulim AND (x + wa1 GT ulim), uct)
  if uct GT 0 then $
    alpha = min([alpha, (ulim(whu)-x(whu))/wa1(whu)])

  ;; Adjust the step according to any boundaries
  wa1 = wa1 * alpha
  wa2 = x + wa1

  ;; If the step put us exactly on a boundary, make sure it is exact
  if lct GT 0 then wa2(whl) = llim(whl)
  if uct GT 0 then wa2(whu) = ulim(whu)

  wa3 = diag * wa1
  pnorm = mpfit_enorm(wa3)

  ;; On the first iteration, adjust the initial step bound
  if iter EQ 1 then delta = min([delta,pnorm])
  
  ;; Evaluate the function at x+p and calculate its norm
  !err = 1
  xnew(ifree) = wa2
  wa4 = call_function(fcn, xnew, _EXTRA=fcnargs)
  nfev = nfev +1
  iflag = !err
  !err = 0
  if iflag LT 0 then begin
      goto, TERMINATE
      errmsg = 'WARNING: premature termination by "'+fcn+'"'
  endif
  fnorm1 = mpfit_enorm(wa4)
  
  ;; Compute the scaled actual reduction
  actred = x(0)*0 - 1.
  if 0.1D * fnorm1 LT fnorm then actred = - (fnorm1/fnorm)^2 + 1.

  ;; Compute the scaled predicted reduction and the scaled directional
  ;; derivative
  for j = 0, n-1 do begin
      wa3(j) = 0;
      wa3(0:j) = wa3(0:j) + fjac(0:j,j)*wa1(ipvt(j))
  endfor

  ;; Remember, alpha is the fraction of the full LM step actually taken
  temp1 = mpfit_enorm(alpha*wa3)/fnorm
  temp2 = (sqrt(alpha*par)*pnorm)/fnorm
  half = temp1(0) * 0 + 0.5
  prered = temp1*temp1 + (temp2*temp2)/half
  dirder = -(temp1*temp1 + temp2*temp2)

  ;; Compute the ratio of the actual to the predicted reduction.
  ratio = half * 0
  tenth = half * 0 + 0.1
  if prered NE 0 then ratio = actred/prered

  ;; Update the step bound
  if ratio LE 0.25D then begin
      if actred GE 0 then temp = half $
      else temp = half*dirder/(dirder + half*actred)
      if ((0.1D*fnorm1) GE fnorm) OR (temp LT 0.1D) then temp = tenth
      delta = temp*min([delta,pnorm/tenth])
      par = par/temp
  endif else begin
      if (par EQ 0) OR (ratio GE 0.75) then begin
          delta = pnorm/half
          par = half*par
      endif
  endelse

  ;; Test for successful iteration
  if ratio GE 0.0001 then begin
      ;; Successful iteration.  Update x, fvec, and their norms
      x = wa2
      wa2 = diag * x

      fvec = wa4
      xnorm = mpfit_enorm(wa2)
      fnorm = fnorm1
      iter = iter + 1
  endif

  ;; Tests for convergence
  if (abs(actred) LE ftol) AND (prered LE ftol) $
    AND  (0.5D * ratio LE 1) then info = 1
  if delta LE xtol*xnorm then info = 2
  if (abs(actred) LE ftol) AND (prered LE ftol) $
    AND (0.5D * ratio LE 1) AND (info EQ 2) then info = 3
  if info NE 0 then goto, TERMINATE

  ;; Tests for termination and stringent tolerances
  if iter GE maxiter then info = 5
  if (abs(actred) LE MACHEP) AND (prered LE MACHEP) $
    AND (0.5*ratio LE 1) then info = 6
  if delta LE MACHEP*xnorm then info = 7
  if gnorm LE MACHEP then info = 8
  if info NE 0 then goto, TERMINATE

  ;; End of inner loop. Repeat if iteration unsuccessful
  if ratio LT 0.0001 then begin
      goto, INNER_LOOP
  endif

  ;; End of outer loop.
  goto, OUTER_LOOP

TERMINATE:
  ;; Termination, either normal or user imposed.
  if iflag LT 0 then info = iflag
  iflag = 0
  if n_elements(ifree) EQ 0 then xnew = xall else xnew(ifree) = x
  if nprint GT 0 then begin
      fvec = call_function(fcn, xnew, _EXTRA=fcnargs)
      fnorm = mpfit_enorm(fvec)
  endif

  fnorm = max([fnorm, fnorm1])
  fnorm = fnorm^2.

  covar = !values.d_nan
  ;; (very carefully) set the covariance matrix COVAR
  if info GT 0 AND NOT keyword_set(nocovar) $
    AND n_elements(n) GT 0 AND n_elements(fvec) GT 0 $
    AND n_elements(fjac) GT 0 AND n_elements(ipvt) GT 0 then begin
      sz = size(fjac)
      if n GT 0 AND sz(0) GT 1 AND sz(1) GE n AND sz(2) GE n $
        AND n_elements(ipvt) GE n then begin
          cv = mpfit_covar(fjac(0:n-1,0:n-1), ipvt(0:n-1))
          nn = n_elements(xall)
          
          ;; Fill in actual covariance matrix, accounting for fixed
          ;; parameters.
          covar = replicate(cv(0)*0, nn, nn)
          for i = 0L, n-1 do begin
              covar(ifree, ifree(i)) = cv(*,i)
          end
          
          ;; Compute errors in parameters
          i = lindgen(nn)
          perror = replicate(covar(0), nn)*0
          wh = where(covar(i,i) GE 0, ct)
          if ct GT 0 then $
            perror(wh) = sqrt(covar(wh, wh))
      endif
  endif
      
  return, xnew
end
FUNCTION Welch,n
;Procedure for calculating Welch Window function (Num Recipes, p 422)
j = findgen(n)
nn = (n-1.)/2.
nd = (n+1.)/2.
RETURN, 1. - ((j - nn)/nd)^2
END ; End Welch 
PRO functc,x,a,f,pder
; Return a function value in f.
; based on lmfunct
; Function for fitting modified Desaubies spectrum
; 23 June 1998 RAV
; F = Ao mu/(1+mu^A2)
; where mu = x/A1

;
;a(0)=F_0
;a(1)=m_*
;a(2)=t+1
;
   mu = x/a(1)				;m/m*
   d = mu^a(2)				;(m/m*)^t
   b = 1./(1+d)
   c = mu*b
   f = a(0)*c				;f0 (m/m*) / (1+(m/m*)^(t+1))
   IF (n_params() GE 4) THEN $
    pder = [ [f/a(0)],  $
             [-f*b*(1+(1-a(2))*d)/a(1)],  $
             [-f*b*alog(mu)*d] ]
   return
END 

FUNCTION mp_funct, p,X=x,Y=y,ERR=sigma
; this function provided to be called by the mpfit.pro
; routine . Devised Aug 12 1998 AC Beresford , Atmos Phys Univ of Adelaide
   nofx = n_elements(x)
   model = fltarr(nofx)
   deriv =fltarr(nofx)
   functc, x,p,model
   deriv = (y-model)/sigma
   return, deriv
END
;
;for one sounding
;
PRO stokes_new,u,v,u_m,v_m,theta,n_bar,f_w,df,z30,no_correct=no_correct
;*****************************************************
;Procedure to calculate Stokes parameters for vertical
;wavenumber spectra
;*****************************************************
;Fourier transform velocity profiles, assumed to be at 
;30 m spacing
; theta is azimuth computed from <u'T90>
   a = size(u)
   n = a(1)
   na = n/2

   is = fltarr(na)
   ds = is & ps = is & qs = is & df = is & phi = is
   df = 0.
   xd = df & yd = df & phi = df & axr = df
      ut = fft(u, -1)
      vt = fft(v, -1)

;Form Stokes parameters following Eckermann and Vincent (1989)
      it = float(ut*conj(ut) + vt*conj(vt))
      dt = float(ut*conj(ut) - vt*conj(vt))
      pt = 2.0*float(conj(ut)*vt)
      qt = 2.0*imaginary(conj(ut)*vt)

      is = it(1:na) & ds = dt(1:na) 
      ps = pt(1:na) & qs = qt(1:na)

      df  = sqrt(total(qt(0:na))^2 + total(pt(0:na))^2 +  $
                    total(dt(0:na))^2)/total(it(0:na))
      phi = atan(total(pt), total(dt))/2.

      xi = asin(abs(total(qt(0:na))/(df*total(it(0:na)))))/2.0
      axr = 1.0/tan(xi)
; Correct for tranverse shear term
   u_t = fltarr(a(1))
   du_t_dz = 0.
   u_t = u_m*cos(theta)-v_m*sin(theta)
;stop
   du_t_dz = total(deriv(z30,u_t))/float(a(1))
;  axr_corr = du_t_dz/n_bar
   axr_corr = du_t_dz/n_bar
   if keyword_set(no_correct) then f_w = axr $
   else f_w = axr-axr_corr
   return
END ; End Stokes

PRO phase_speed_new,u_m,v_m,phi,up,vp,f,n_bar,f_w,m_bar,k_bar, $
                omega,c_z,c_i,c_x,c_y,df,c_ix,c_iy,c_gx,c_gy, $
;               c_xcomp,c_ycomp
                c_xcomp,c_ycomp, u_mean, u_mm, v_mm, theta,dzm
;stop
; Procedure to compute phase speeds
   a = size(up)
;  b = where(f_w GT 1 AND f_w LE 10, n_good)
;compute azimuth of mean winda = size(up)
   u_mm = reform(rebin(u_m,1))
   v_mm = reform(rebin(v_m,1))
   theta = atan(u_mm,v_mm)
   u_mean = sqrt(u_mm^2+v_mm^2)
;  mean_spec,up,vp,m_bar,m,lz,1,30.
;stop
   mean_spec_new,up,vp,m_bar,m,lz,1,dzm

   m_z = 2.0*!pi*m_bar

; First compute horizontal wavenumber
   k_bar = 0.
   k_bar = f*m_z*sqrt(f_w^2-1)/n_bar
   l_x = 2.0*!pi/k_bar

; Find intrinsic phase speeds
   c_z = f*f_w/(m_z)
   c_i = f*f_w/k_bar
   c_ix = c_i*sin(phi)
   c_iy = c_i*cos(phi)

; Now compute group velocity c_g = c_i + U
   c_gx = c_ix + u_mm
   c_gy = c_iy + v_mm
;plot_cg,c_gx,c_gy

; Find ground based phase speed, c
; and ground-based frequency
   omega = f*f_w + k_bar*u_mean*cos(theta-phi)
   c =  abs((c_i + u_mean*cos(theta-phi)))
   c_x = c/sin(phi)      ; note division since not true component velocity
   c_y = c/cos(phi)
   c_xcomp = (c_i + u_mean*cos(theta-phi))*sin(phi)
   c_ycomp = (c_i + u_mean*cos(theta-phi))*cos(phi)

   return
END                             ; End of phase speed

;using a complete dispersion relation
PRO phase_speed_new2,u_m,v_m,phi,up,vp,f,n_bar,f_w,m_bar,k_bar, $
                omega,c_z,c_i,c_x,c_y,df,c_ix,c_iy,c_gx,c_gy, $
;               c_xcomp,c_ycomp
                c_xcomp,c_ycomp, u_mean, u_mm, v_mm, theta,dzm,alpha
;stop
; Procedure to compute phase speeds
   a = size(up)
;  b = where(f_w GT 1 AND f_w LE 10, n_good)
;compute azimuth of mean winda = size(up)
   u_mm = reform(rebin(u_m,1))
   v_mm = reform(rebin(v_m,1))
   theta = atan(u_mm,v_mm)
   u_mean = sqrt(u_mm^2+v_mm^2)
;  mean_spec,up,vp,m_bar,m,lz,1,30.
;stop
   mean_spec_new,up,vp,m_bar,m,lz,1,dzm

   m_z = 2.0*!pi*m_bar

; First compute horizontal wavenumber
;  k_bar = f*m_z*sqrt(f_w^2-1)/n_bar
;  if abs(f*f_w) ge n_bar then message,'frequency too high'
   k_bar = f*sqrt(m_z^2+alpha^2)*sqrt(f_w^2-1)/sqrt(n_bar^2-(f*f_w)^2)
   l_x = 2.0*!pi/k_bar
;print,sqrt(f_w^2-1)
;print,sqrt(n_bar^2-(f*f_w)^2)
;print,k_bar
;print,k_bar
;print,l_x/1000.

; Find intrinsic phase speeds
   c_z = f*f_w/(m_z)
   c_i = f*f_w/k_bar
   c_ix = c_i*sin(phi)
   c_iy = c_i*cos(phi)

; Now compute group velocity c_g = c_i + U
   c_gx = c_ix + u_mm
   c_gy = c_iy + v_mm
;plot_cg,c_gx,c_gy

; Find ground based phase speed, c
; and ground-based frequency
   omega = f*f_w + k_bar*u_mean*cos(theta-phi)
   c =  abs((c_i + u_mean*cos(theta-phi)))
   c_x = c/sin(phi)      ; note division since not true component velocity
   c_y = c/cos(phi)
   c_xcomp = (c_i + u_mean*cos(theta-phi))*sin(phi)
   c_ycomp = (c_i + u_mean*cos(theta-phi))*cos(phi)

   return
END                             ; End of phase speed

;************************************************************
;
;************************************************************
;
;PRO directions,up,vp,tp,rho,phi,ut,vt,uq,ui
;PRO directions,up,vp,tp,rho,phi,ut,vt,uq,vq,ui,vi	;newly modified
PRO directions_new,up,vp,tp,rho,phi,ut,vt,uq,vq,ui,vi	;newly modified
; Procedure to compute mean direction of wave propagation
; from velocity and normalized temperature perturbations.
; Use covariances <u'T'+90> and <v'T'+90> to determine
; azimuth from North. Where T'+90 is the value of
; temperature after shifting phase by +90 deg via Hilbert
; transform.

   a=size(up)
   phi = 0.					;dir of wave propagation
   ut = phi & vt = phi
   uq = ut & ui = ut
   vq = ut						;newly added
; Carry out Hilbert transform
      to = float(HILBERT(tp, -1))

      ti = tp
      x = total(to*up*rho)/float(a(1))	;column average
      y = total(to*vp*rho)/float(a(1))	;column average
      ut = x
      vt = y
      uq = total(to*up)/n_elements(to)
      vq = total(to*vp)/n_elements(to)		;newly added
      ui = total(ti*up)/n_elements(ti)
; Phi is azimuth from north
      phi = atan(x, y)
   return
END                             ; End Directions

;*************************************************
PRO t_spec_new,t,tspec,err,m,lamda,navr,cpm
; Procedure to calulate vertical-wavenumber power spectrum for temperature
; m and lamda are the output vertical wavenumber and scale arrays
; averaged over navr points
; cpm is the sample spacing in metres
; tspec is the output array (real) of length npts/navr
; err is a standard deviation of the tspec at each wavenumber
; Generate window function
; Average over all available observations if necessary
   a=size(t)
   npts=a(1)
;  IF (a(0) EQ 2) THEN n_obs = a(2) $
;   ELSE n_obs = 1
   n_obs = 1
   window = welch(npts)
   nred = npts/navr
   ts = fltarr(npts)
   tsq =fltarr(npts)
   nmax=npts-1
   nf = float(npts)/(total(window*window))

   delm = 1.0/(cpm*float(nred))
   m = findgen(nred)*delm
   if (navr ne 1) then m = m + delm/2.		; ? my question
   m(0) = delm/2.
   uf = complex((t(0:nmax)-total(t(0:nmax))/npts)*window ,0.0)
   ut = fft(uf,-1)
   spec = float(ut*conj(ut))
   ts = ts +spec
   tsq =spec*spec+tsq
;stop,'777'

   fnp =float(n_obs)
   norm = float(2*navr)*nf/(n_obs*delm) ; NB One-sided spectrum
   err = tsq/fnp - ts*ts/(fnp*fnp)
   err = fnp*sqrt(err)
   tspec = rebin(ts,nred)*norm
   tspec = tspec(0:npts/2)
;  err = rebin(err,nred)*norm
   err = 0.1*abs(ts)*norm			;I add this
   err = err(0:npts/2)
   lamda = 1.0/(m*1000)

return
END ; t_spec

PRO mean_spec_new,u,v,m_bar,m,lamda,navr,cpm
;Procedure to calulate mean vertical-wavenumber for wind components
;m and lamda are the output vertical wavenumber and wavelength arrays
;averaged over navr points
;cpm is the sample spacing in metres
;uspec and vspec are the output arrays (real) of length npts/navr
;Generate window function
;m_bar is the number of mean wavenumbers
   a=size(u)
   npts=a(1)
   n2 = npts/2
;  IF (a(0) EQ 2) THEN n_data = a(2)
;  m_bar = reform(fltarr(n_data))
   m_bar = 0.
   window = welch(npts)
   nred = npts/navr
   nmax=npts-1
   nf = float(npts)/(total(window*window))

   delm = 1.0/(cpm*float(nred))
   m = findgen(nred)*delm
   IF (navr NE 1) THEN  m = m + delm/2.
   m(0) = delm/2. 
   norm = float(navr)*nf/delm

   i = complex(0,1)
;  IF (n_data GT 1) THEN BEGIN
;     FOR i = 0, n_data-1 DO BEGIN
         uf = (u(0:nmax)-total(u(0:nmax))/npts)*window
         vf = (v(0:nmax)-total(v(0:nmax))/npts)*window
         ut = fft(uf, -1)
         us = float(ut*conj(ut))
         vt = fft(vf, -1)
         vs = float(vt*conj(vt))
         us =  us + vs

         if (navr ne 1) THEN  $
          us = rebin(us, nred)*norm  $
         ELSE  $
          us = us*norm

         m_bar = rebin(us(1:n2)*m(1:n2),1)/rebin(us(1:n2),1)
;     ENDFOR
;  ENDIF 
   lamda  = 1.0/(m*1000.0)
   RETURN
END
;*********************************************************
PRO tspec_plot_new, ts, m,lamda,ttl=ttl
   a = size(ts)
   n = a(1)/2-10
;   plot_oo,m(1:n),ts(1:n),xstyle = 9, $
;    xmargin =[15,0],$
;    xtitle='!17WAVENUMBER (cpm)',$
;;   ytitle='!17Normalized T PSD (cpm!E-1!N)',/NODATA
;    ytitle='!17Norma T PSD (cpm!E-1!N)',/NODATA
;   oplot,m(1:n),ts(1:n)
;   if keyword_set(ttl) then xttl=ttl+'!C'+'WAVELENGTH (km)' $
;   else xttl='WAVELENGTH (km)'
;   axis,xaxis=1,xstyle = 1, $
;    xrange=[lamda(1), $
;            lamda(n)],xtitle=xttl
   plot_oo,m(1:n),ts(1:n),xstyle = 9, $
    xmargin =[15,0],/NODATA
   oplot,m(1:n),ts(1:n)
   axis,xaxis=1,xstyle = 1, $
    xrange=[lamda(1), $
            lamda(n)],xtitle=ttl
   return
END

PRO directions_new_bandpass,up,vp,tp,rho,phi
; Procedure to compute mean direction of wave propagation
; from velocity and normalized temperature perturbations.
; Use covariances <u'T'+90> and <v'T'+90> to determine
; azimuth from North. Where T'+90 is the value of
; temperature after shifting phase by +90 deg via Hilbert
; transform.

   a=size(up)
   phi = 0.					;dir of wave propagation
   ut = phi & vt = phi
   uq = ut & ui = ut
   vq = ut						;newly added
; Carry out Hilbert transform
      to = float(HILBERT(tp, -1))

      ti = tp
      x = total(to*up*rho)/float(a(1))	;column average
      y = total(to*vp*rho)/float(a(1))	;column average
; Phi is azimuth from north
      phi = atan(x, y)
   return
END                             ; End Directions

pro get_perturbation,up,vp,tp,tb,um,vm,npoly
s=size(up)

zkm=findgen(s(1))*0.03

;--> uraw,vraw,traw
for i=0,s(2)-1 do begin
    if up(0,i) ne -999. and up(0,i) ne 999. then begin
        uraw=reform(up(*,i)+um(*,i))
        vraw=reform(vp(*,i)+vm(*,i))
        co = poly_fit(double(zkm),double(uraw),npoly,fitu)
        co = poly_fit(double(zkm),double(vraw),npoly,fitv)
        um(*,i)=fitu
        vm(*,i)=fitv
        up(*,i)=uraw-fitu
        vp(*,i)=vraw-fitv
    endif
    if tp(0,i) ne -999. and tp(0,i) ne 999. then begin
        traw=reform(tp(*,i)+tb(*,i))
        co = poly_fit(double(zkm),double(traw),npoly,fitt)
        tb(*,i)=fitt
        tp(*,i)=traw-fitt
    endif
endfor

return
end

pro get_perturbation_ltm,up,vp,tp,tb,um,vm,mn,ultm,vltm,tltm,npoly
s=size(up)

zkm=findgen(s(1))*0.03

;--> uraw,vraw,traw
for i=0,s(2)-1 do begin
    j=mn(i)-1
    if up(0,i) ne -999. and up(0,i) ne 999. then begin
        if ultm(0,j) eq -999. then stop
        if vltm(0,j) eq -999. then stop
        uraw=reform(up(*,i)+um(*,i))-ultm(*,j)
        vraw=reform(vp(*,i)+vm(*,i))-vltm(*,j)
        co = poly_fit(double(zkm),double(uraw),npoly,fitu)
        co = poly_fit(double(zkm),double(vraw),npoly,fitv)
        um(*,i)=fitu+ultm(*,j)
        vm(*,i)=fitv+vltm(*,j)
        up(*,i)=uraw-fitu
        vp(*,i)=vraw-fitv
    endif
    if tp(0,i) ne -999. and tp(0,i) ne 999. then begin
        if tltm(0,j) eq -999. then stop
        traw=reform(tp(*,i)+tb(*,i))-tltm(*,j)
        co = poly_fit(double(zkm),double(traw),npoly,fitt)
        tb(*,i)=fitt+tltm(*,j)
        tp(*,i)=traw-fitt
    endif
endfor

return
end

pro get_high_pass,up,vp,tp,tb,um,vm,bndy
s=size(up)

zkm=findgen(s(1))*0.03

fhigh=1
flow=2*0.03/bndy

;--> uraw,vraw,traw
for i=0,s(2)-1 do begin
    if up(0,i) ne -999. and up(0,i) ne 999. then begin
        uraw=reform(up(*,i)+um(*,i))
        vraw=reform(vp(*,i)+vm(*,i))
        co = poly_fit(double(zkm),double(uraw),1,fitu)
        co = poly_fit(double(zkm),double(vraw),1,fitv)
        uraw=uraw-fitu
        vraw=vraw-fitv
        um(*,i)=digital_smooth(uraw,flow,fhigh)+fitu
        vm(*,i)=digital_smooth(vraw,flow,fhigh)+fitv
        up(*,i)=uraw-um(*,i)+fitu
        vp(*,i)=vraw-vm(*,i)+fitv
    endif
    if tp(0,i) ne -999. and tp(0,i) ne 999. then begin
        traw=reform(tp(*,i)+tb(*,i))
        co = poly_fit(double(zkm),double(traw),1,fitt)
        traw=traw-fitt
        tb(*,i)=digital_smooth(traw,flow,fhigh)+fitt
        tp(*,i)=traw-tb(*,i)+fitt
    endif
endfor

return
end

;modified on 2/27/04 to add the option to change all the alog() to alog10()
;i.e., add fit2_new, fit4_new, and fit5_new

;*********************************************************
;
;follows what's used in Nastrom and VanZandt, 2001: spectral slopes at
;high wavenumber will be determined over the interval 6/6400 =< m <= 27/6400
;see the function 'fit2'
;
pro fit2,ts,m,t
m0=6./6400.
m1=27./6400.
ava=where(m ge m0 and m le m1)
logts=alog(ts(ava))
logm=alog(m(ava))
result=poly_fit(logm,logts,1,/double)
t=-result(1)

return
end

pro fit2_new,ts,m,t
m0=6./6400.
m1=27./6400.
ava=where(m ge m0 and m le m1)
logts=alog10(ts(ava))
logm=alog10(m(ava))
result=poly_fit(logm,logts,1,/double)
t=-result(1)

return
end

PRO fit3, ts,sigma,m,t,a
;
;modified from fit
;exactly the same as fit but letting the initial guess of t being t2 (derived from fit2)
;
; Procedure to fit m-spectrum for normalized temperature spectrum (ts)
;; to modified Desaubies spectrum in vertical wavenumber (m)

   dummystring = ' '
   npt = n_elements(m)          ; Make spectrum one-sided
   n = npt/2.-1
   s = ts(1:n)
   x = m(1:n)
   ss = sigma(1:n)
; Make initial guess of fitting constants
; Get mean-square normalized temperature
   dm = m(4)-m(3)
   ts_bar = total(s)*dm
; Next mean wavenumber
   m_bar = total(s*x)*dm/ts_bar
; First guess at mr/ms
;  t = 4.0
  a = [fo,ms,t ]
   ms = x(1)
   fo = ts_bar/ms
   a = [fo,ms,t+1 ]
   fa = {X:x, Y:s, ERR:ss}
   functc,x,a,f        ; f is starting model
   deviates = (s-f)
;fitting a model to the avaerage spectrum using fitting algorithm
;embodied in mpfit.pro
   p_fit = mpfit('mp_funct',a,functargs= fa,/QUIET)
   functc,x,p_fit,spec_fit      ; evaluate fitted function
;loadct,34
   Oplot,x,spec_fit,line = 2,color=74    ; plot it over data
loadct,0
   a = p_fit
   return
END

;
;same as fit2, but 1./1500. =< m <= 1./300.
;
pro fit4,ts,m,t
m0=1./1500.
m1=1./300.
ava=where(m ge m0 and m le m1)
logts=alog(ts(ava))
logm=alog(m(ava))
result=poly_fit(logm,logts,1,/double)
t=-result(1)

return
end

pro fit4_new,ts,m,t
m0=1./1500.
m1=1./300.
ava=where(m ge m0 and m le m1)
logts=alog10(ts(ava))
logm=alog10(m(ava))
result=poly_fit(logm,logts,1,/double)
t=-result(1)

return
end

;same as fit2, but 1./1000. =< m <= 1./200.
pro fit5,ts,m,t
m0=1./1000.
m1=1./200.
ava=where(m ge m0 and m le m1)
logts=alog(ts(ava))
logm=alog(m(ava))
result=poly_fit(logm,logts,1,/double)
t=-result(1)

return
end

pro fit5_new,ts,m,t
m0=1./1000.
m1=1./200.
ava=where(m ge m0 and m le m1)
logts=alog10(ts(ava))
logm=alog10(m(ava))
result=poly_fit(logm,logts,1,/double)
t=-result(1)

return
end

PRO fit6, ts,sigma,m,a
; exactly the same as fit, but with t already being estimated previously
; and only fit ms and f0

   dummystring = ' '
   npt = n_elements(m)          ; Make spectrum one-sided
   n = npt/2.-1
   s = ts(1:n)
   x = m(1:n)
   ss = sigma(1:n)
; Make initial guess of fitting constants
; Get mean-square normalized temperature
   dm = m(4)-m(3)
   ts_bar = total(s)*dm
; Next mean wavenumber
   m_bar = total(s*x)*dm/ts_bar
; First guess at m_star
   ms = x(1)
   fo = ts_bar/ms

;  t = 4.0
;  a = [fo,ms,t ]
   a = [fo,ms]

   fa = {X:x, Y:s, ERR:ss}

;  functc,x,a,f        ; f is starting model
   functc2,x,a,f        ; f is starting model

   deviates = (s-f)
; fitting a model to the avaerage spectrum using fitting algorithm
; embodied in mpfit.pro

   p_fit = mpfit('mp_funct2',a,functargs= fa,/QUIET)

;  functc,x,p_fit,spec_fit      ; evaluate fitted function
   functc2,x,p_fit,spec_fit      ; evaluate fitted function

;loadct,6
;   Oplot,x,spec_fit,line = 2,color=210    ; plot it over data
;loadct,0
   a = p_fit
   return
END

PRO fit10, ts,m,ms,t
; fit ts to the form F(m) = F0 (m/m*) / (1+(m/m*)^(t+1))
; using curvefit.pro from IDL

   itmax=200
   dummystring = ' '
   npt = n_elements(m)          ; Make spectrum one-sided
   n = npt/2.-1
   s = ts(1:n)
   x = m(1:n)
; Make initial guess of fitting constants
; Get mean-square normalized temperature
   dm = m(4)-m(3)
   ts_bar = total(s)*dm
; Next mean wavenumber
   m_bar = total(s*x)*dm/ts_bar
; First guess at m_star
   ms = x(1)
   fo = ts_bar/ms

   t0 = 4.0
   a=[fo,ms,t0]
   weights=1./s
   result=curvefit(x,s,weights,a,iter=iter,function_name='functc',itmax=itmax,tol=1.0d-4)
   ms=a(1)
   t=a(2)-1.
   if iter eq itmax then begin
       ms=-999.
       t=-999.
   endif
   return
END

PRO fit11, ts,m,ms
; same as fit10 but with t already known from fit5 (t5)
; using curvefit.pro from IDL

   itmax=200
   dummystring = ' '
   npt = n_elements(m)          ; Make spectrum one-sided
   n = npt/2.-1
   s = ts(1:n)
   x = m(1:n)
; Make initial guess of fitting constants
; Get mean-square normalized temperature
   dm = m(4)-m(3)
   ts_bar = total(s)*dm
; Next mean wavenumber
   m_bar = total(s*x)*dm/ts_bar
; First guess at m_star
   ms = x(1)
   fo = ts_bar/ms

   a=[fo,ms]
   weights=1./s
   result=curvefit(x,s,weights,a,iter=iter,function_name='functc2',itmax=itmax,tol=1.0d-4)
   ms=a(1)
;print,iter
;  if iter eq itmax then ms=-999.
   return
END

PRO mean_spec2,u,v,m_bar,m,lamda,navr,cpm
;
; exactly the same as mean_spec but calling welch2 instead of welch
;
;Procedure to calulate mean vertical-wavenumber for wind components
;m and lamda are the output vertical wavenumber and wavelength arrays
;averaged over navr points
;cpm is the sample spacing in metres
;uspec and vspec are the output arrays (real) of length npts/navr
;Generate window function
;m_bar is the number of mean wavenumbers
   a=size(u)
   npts=a(1)
   n2 = npts/2
   IF (a(0) EQ 2) THEN n_data = a(2)
   m_bar = reform(fltarr(n_data))
   window = welch2(npts)
   nred = npts/navr
   nmax=npts-1
   nf = float(npts)/(total(window*window))

   delm = 1.0/(cpm*float(nred))
   m = findgen(nred)*delm
   IF (navr NE 1) THEN  m = m + delm/2.
   m(0) = delm/2. 
   norm = float(navr)*nf/delm

   i = complex(0,1)
   IF (n_data GT 1) THEN BEGIN
      FOR i = 0, n_data-1 DO BEGIN
         uf = (u(0:nmax,i)-total(u(0:nmax,i))/npts)*window
         vf = (v(0:nmax,i)-total(v(0:nmax,i))/npts)*window
         ut = fft(uf, -1)
         us = float(ut*conj(ut))
         vt = fft(vf, -1)
         vs = float(vt*conj(vt))
         us =  us + vs

         if (navr ne 1) THEN  $
          us = rebin(us, nred)*norm  $
         ELSE  $
          us = us*norm

         m_bar(i) = rebin(us(1:n2)*m(1:n2),1)/rebin(us(1:n2),1)
      ENDFOR
   ENDIF 
   lamda  = 1.0/(m*1000.0)
   RETURN
END

PRO t_spec2,t,tspec,err,m,lamda,navr,cpm
;
; exactly the same as t_spec but calling welch2 instead of welch
;
; Procedure to calulate vertical-wavenumber power spectrum for temperature
; m and lamda are the output vertical wavenumber and scale arrays
; averaged over navr points
; cpm is the sample spacing in metres
; tspec is the output array (real) of length npts/navr
; err is a standard deviation of the tspec at each wavenumber
; Generate window function
; Average over all available observations if necessary
   a=size(t)
   npts=a(1)
   IF (a(0) EQ 2) THEN n_obs = a(2) $
    ELSE n_obs = 1
   window = welch2(npts)
   nred = npts/navr
   ts = fltarr(npts)
   tsq =fltarr(npts)
   nmax=npts-1
   nf = float(npts)/(total(window*window))

   delm = 1.0/(cpm*float(nred))
   m = findgen(nred)*delm
   if (navr ne 1) then m = m + delm/2.		; ? my question
   m(0) = delm/2.
   FOR i = 0,n_obs-1 DO BEGIN
      uf = complex((t(0:nmax,i)-total(t(0:nmax,i))/npts)*window ,0.0)
      ut = fft(uf,-1)
      spec = float(ut*conj(ut))
      ts = ts +spec
      tsq =spec*spec+tsq
   ENDFOR
   fnp =float(n_obs)
   norm = float(2*navr)*nf/(n_obs*delm) ; NB One-sided spectrum
   err = tsq/fnp - ts*ts/(fnp*fnp)
   err = fnp*sqrt(err)
   tspec = rebin(ts,nred)*norm
   tspec = tspec(0:npts/2)
   err = rebin(err,nred)*norm
   err = err(0:npts/2)
   lamda = 1.0/(m*1000)
return
END ; t_spec

PRO rotary_spec2,u,v,acp,ccp,m,lamda,navr,cpm
;
; exactly the same as rotary_spec but call welch2 instead of welch
;
;Procedure to calulate vertical-wavenumber power spectrum for wind components
;m and lamda are the output vertical wavenumber and wavelength arrays
;averaged over navr points
;cpm is the sample spacing in metres
;uspec and vspec are the output arrays (real) of length npts/navr

;Generate window function
   a=size(u)
   npts=a(1)
   n2 = npts/2
   IF (a(0) EQ 2) THEN n_data = a(2) ELSE n_DATA =1
   window = welch2(npts)
   nred = npts/navr
   nmax=npts-1
   nf = float(npts)/(total(window*window))
   delm = 1.0/(cpm*float(nred))
   m = findgen(nred)*delm
   IF (navr NE 1) THEN m = m + delm/2.
   m(0) = delm/2.
   norm = float(navr)*nf/delm
   acp = fltarr(n2) & ccp = acp

   FOR i = 0,n_data-1 DO BEGIN
      uf = (u(0:nmax,i)-total(u(0:nmax,i))/npts)*window
      vf = (v(0:nmax,i)-total(v(0:nmax,i))/npts)*window
      ut = fft(complex(uf,vf), -1)
      us = float(ut*conj(ut))
      IF (navr NE 1) THEN  $
       us = rebin(us,nred)*norm ELSE us = us*norm
;Anticlockwise component is positive half space
;clockwise component is negative half space
      acp = us(1:n2)+acp
      ccp = rotate(us(n2:nmax),2)+ccp
   ENDFOR
   acp = acp/float(n_data)
   ccp = ccp/float(n_data)
   lamda  = 1.0/(m*1000.0)
   return
END                             ; End rotary_spec
PRO phase_speed2,u_m,v_m,phi,up,vp,f,n_bar,f_w,m_bar,k_bar, $
                omega,c_z,c_i,c_x,c_y,df,c_ix,c_iy,c_gx,c_gy, $
;               c_xcomp,c_ycomp
                c_xcomp,c_ycomp, u_mean, u_mm, v_mm, theta
;
; exactly the same as phase_speed but call mean_spec2 instead of mean_spec
;
; Procedure to compute phase speeds
   a = size(up)
   b = where(f_w GT 1 AND f_w LE 10, n_good)
;compute azimuth of mean winda = size(up)
   u_mm = reform(rebin(u_m,1,a(2)))
   v_mm = reform(rebin(v_m,1,a(2)))
   theta = atan(u_mm,v_mm)
   u_mean = sqrt(u_mm^2+v_mm^2)
   mean_spec2,up,vp,m_bar,m,lz,1,30.

   m_z = 2.0*!pi*m_bar

; First compute horizontal wavenumber
   k_bar = fltarr(n_good)
   k_bar = f*m_z(b)*sqrt(f_w(b)^2-1)/n_bar
   l_x = 2.0*!pi/k_bar

; Find intrinsic phase speeds
   c_z = f*f_w(b)/(m_z(b))
   c_i = f*f_w(b)/k_bar
   c_ix = c_i*sin(phi(b))
   c_iy = c_i*cos(phi(b))

; Now compute group velocity c_g = c_i + U
   c_gx = c_ix + u_mm(b)
   c_gy = c_iy + v_mm(b)
;plot_cg,c_gx,c_gy

; Find ground based phase speed, c
; and ground-based frequency
   omega = f*f_w(b) + k_bar*u_mean(b)*cos(theta(b)-phi(b))
   c =  abs((c_i + u_mean(b)*cos(theta(b)-phi(b))))
   c_x = c/sin(phi(b))      ; note division since not true component velocity
   c_y = c/cos(phi(b))
   c_xcomp = (c_i + u_mean(b)*cos(theta(b)-phi(b)))*sin(phi(b))
   c_ycomp = (c_i + u_mean(b)*cos(theta(b)-phi(b)))*cos(phi(b))

   df = df(b)

   return
END                             ; End of phase speed

; postdarkening (Nastrom and VanZandt (2001) p. 14,370)
pro postdark, tsh, m
n=n_elements(m)
j=1.+findgen(n)
tsh=tsh/(2.0*(1-cos(2.0*!pi*j/n)))
return
end
;****************************************************************
FUNCTION pressure,z
; Procedure to produce pressure  profile
; for height range from 18 - 25 km
   zO = 0.250
;Standard Height-Pressure levels for Radiosonde Data at Cocos Island
   hstd = [8.84, 8.77, 8.70, 8.57, 8.32, 8.16,  $
           7.96, 7.83, 7.66, 7.44, $
           7.14, 6.96, 6.85, 6.75, 6.72]
   zstd = [0.100, 0.784, 1.514, 3.157, 5.869,  $
           7.581, 9.686, 10.955, 12.432, $ 
           14.219, 16.553, 18.603, 20.617, 23.766, 26.389]

   pstd = [1000., 925, 850, 700, 500, 400, 300,  $
           250, 200, 150, 100, $
           70, 50, 30, 20]

   p = fltarr(n_elements(z))
   po = interpol(pstd,zstd,zO)
   po = po(0)
   H = interpol(hstd,zstd,z)
   p = po*exp(-(z-zO)/H)
; Convert pressure to Pascals
   p = p*100.
   return,p
END 
PRO mean_spec,u,v,m_bar,m,lamda,navr,cpm
;Procedure to calulate mean vertical-wavenumber for wind components
;m and lamda are the output vertical wavenumber and wavelength arrays
;averaged over navr points
;cpm is the sample spacing in metres
;uspec and vspec are the output arrays (real) of length npts/navr
;Generate window function
;m_bar is the number of mean wavenumbers
   a=size(u)
   npts=a(1)
   n2 = npts/2
   IF (a(0) EQ 2) THEN n_data = a(2)
   m_bar = reform(fltarr(n_data))
   window = welch(npts)
   nred = npts/navr
   nmax=npts-1
   nf = float(npts)/(total(window*window))

   delm = 1.0/(cpm*float(nred))
   m = findgen(nred)*delm
   IF (navr NE 1) THEN  m = m + delm/2.
   m(0) = delm/2. 
   norm = float(navr)*nf/delm

   i = complex(0,1)
   IF (n_data GT 1) THEN BEGIN
      FOR i = 0, n_data-1 DO BEGIN
         uf = (u(0:nmax,i)-total(u(0:nmax,i))/npts)*window
         vf = (v(0:nmax,i)-total(v(0:nmax,i))/npts)*window
         ut = fft(uf, -1)
         us = float(ut*conj(ut))
         vt = fft(vf, -1)
         vs = float(vt*conj(vt))
         us =  us + vs

         if (navr ne 1) THEN  $
          us = rebin(us, nred)*norm  $
         ELSE  $
          us = us*norm

         m_bar(i) = rebin(us(1:n2)*m(1:n2),1)/rebin(us(1:n2),1)
      ENDFOR
   ENDIF 
   lamda  = 1.0/(m*1000.0)
   RETURN
END
;*********************************************************************
PRO stokes,u,v,u_m,v_m,theta,n_bar,f_w,df
;*****************************************************
;Procedure to calculate Stokes parameters for vertical
;wavenumber spectra
;*****************************************************
;Fourier transform velocity profiles, assumed to be at 
;30 m spacing
; theta is azimuth computed from <u'T90>
   a = size(u)
   n = a(1)
   na = n/2
   loop =a(2)
   z30 = findgen(a(1))*30+18000.

   is = fltarr(na, loop)
   ds = is & ps = is & qs = is & df = is & phi = is
   df = fltarr(loop)
   xd = df & yd = df & phi = df & axr = df
   FOR k = 0, loop -1 DO BEGIN
      ut = fft(u(*, k), -1)
      vt = fft(v(*, k), -1)

;Form Stokes parameters following Eckermann and Vincent (1989)

      it = float(ut*conj(ut) + vt*conj(vt))
      dt = float(ut*conj(ut) - vt*conj(vt))
      pt = 2.0*float(conj(ut)*vt)
      qt = 2.0*imaginary(conj(ut)*vt)

      is(*, k) = it(1:na) & ds(*, k) = dt(1:na) 
      ps(*, k) = pt(1:na) & qs(*, k) = qt(1:na)

      df(k)  = sqrt(total(qt(0:na))^2 + total(pt(0:na))^2 +  $
                    total(dt(0:na))^2)/total(it(0:na))
      phi(k) = atan(total(pt), total(dt))/2.

      xi = asin(abs(total(qt(0:na))/(df(k)*total(it(0:na)))))/2.0
      axr(k) = 1.0/tan(xi)
   ENDFOR
; Correct for tranverse shear term
   u_t = fltarr(a(1))
   du_t_dz = fltarr(a(2))
   FOR i = 0,a(2)-1 DO BEGIN 
      u_t = u_m(*,i)*cos(theta(i))-v_m(*,i)*sin(theta(i))
      du_t_dz(i) = total(deriv(z30,u_t))/float(a(1))
   ENDFOR 
   axr_corr = du_t_dz/n_bar
   f_w = axr-axr_corr
END ; End Stokes
;*************************************************************
PRO phase_speed,u_m,v_m,phi,up,vp,f,n_bar,f_w,m_bar,k_bar, $
                omega,c_z,c_i,c_x,c_y,df,c_ix,c_iy,c_gx,c_gy, $
;               c_xcomp,c_ycomp
                c_xcomp,c_ycomp, u_mean, u_mm, v_mm, theta
; Procedure to compute phase speeds
   a = size(up)
   b = where(f_w GT 1 AND f_w LE 10, n_good)
;compute azimuth of mean winda = size(up)
   u_mm = reform(rebin(u_m,1,a(2)))
   v_mm = reform(rebin(v_m,1,a(2)))
   theta = atan(u_mm,v_mm)
   u_mean = sqrt(u_mm^2+v_mm^2)
   mean_spec,up,vp,m_bar,m,lz,1,30.

   m_z = 2.0*!pi*m_bar

; First compute horizontal wavenumber
   k_bar = fltarr(n_good)
   k_bar = f*m_z(b)*sqrt(f_w(b)^2-1)/n_bar
   l_x = 2.0*!pi/k_bar

; Find intrinsic phase speeds
   c_z = f*f_w(b)/(m_z(b))
   c_i = f*f_w(b)/k_bar
   c_ix = c_i*sin(phi(b))
   c_iy = c_i*cos(phi(b))

; Now compute group velocity c_g = c_i + U
   c_gx = c_ix + u_mm(b)
   c_gy = c_iy + v_mm(b)
;plot_cg,c_gx,c_gy

; Find ground based phase speed, c
; and ground-based frequency
   omega = f*f_w(b) + k_bar*u_mean(b)*cos(theta(b)-phi(b))
   c =  abs((c_i + u_mean(b)*cos(theta(b)-phi(b))))
   c_x = c/sin(phi(b))      ; note division since not true component velocity
   c_y = c/cos(phi(b))
   c_xcomp = (c_i + u_mean(b)*cos(theta(b)-phi(b)))*sin(phi(b))
   c_ycomp = (c_i + u_mean(b)*cos(theta(b)-phi(b)))*cos(phi(b))

   df = df(b)

   return
END                             ; End of phase speed
PRO angular_spec,phi,en,title,mean_dir
; Procedure to plot angular spectrum in polar form
; First compute histogram of directions
; Phi assumed to be clockwise from north
; En is the wave energy for each profile
   dir = phi*!radeg
   b = where(dir lt 0,n_dir)
   IF (n_dir GT 0) THEN dir(b)=dir(b)+360.

   r_d = histogram(dir,bin = 30,min = 0,max = 360,reverse_indices = h)
   a = size(r_d)
   r = fltarr(a(1))
   theta = findgen(13)*30*!dtor

; h gives indices of elements in each bin
; e_tot is the total energy for the period under consideration
   e_tot = total(en)
   FOR i = 0,11 DO BEGIN
      IF (r_d(i) NE 0) THEN BEGIN
         v = h(h(i):h(i+1)-1)
         r(i) = total(en(v))/e_tot ;Angular spectrum
      ENDIF ELSE BEGIN
         r(i) = 0
      ENDELSE
   ENDFOR
 goto,jump4
   !x.style = 1
   !y.style = 1
   !x.range = [-0.5,0.5]
   !y.range = [-0.5,0.5]
   !x.tickv = [-0.5,-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3,0.4, 0.5]
   !x.tickname = ['!170.4','0.2', ' ', '0.2', '0.4']
   !y.tickv = [ -0.5 ,-0.4,-0.3,-0.2,-0.1, 0.0, 0.1, 0.2, 0.3,0.4, 0.5]
   !y.tickname = ['0.4','0.2',' ','0.2','0.4']
   df = 30.0*!dtor
;  the array atheta allows the plotting in polar co-ords to
; plot with N to top of page, and go clockwise in angle
   atheta = !pi*.5 - theta
   b = where(atheta lt 0.)
   atheta(b) = atheta(b) + 2.*!pi
   plot,r,atheta,xst = 5,yst = 5,/polar,/nodata,$
    title = title,charsize = 1.5 ,$
    xmargin = [5,5]
   oplot,[0,r(0)],[atheta(0),atheta(0)],/polar
   oplot,[0,r(0)],[atheta(0)-df,atheta(0)-df],/polar
   arc_circ,r(0),atheta(0),-df
   FOR i = 1,11 DO BEGIN
      oplot,[0,r(i)],[atheta(i),atheta(i)],/polar
      oplot,[0,r(i)],[atheta(i)-df,atheta(i)-df],/polar
      arc_circ,r(i),atheta(i),-df
   ENDFOR
   axis,0,0,xax = 0,charsize = 1.3
   axis,0,0,yax = 0,charsize = 1.3
   xyouts,0.495,0.02,'E',charsize = 1.0
   xyouts,0.015,0.45,'N',charsize = 1.0
; Reset plot parameters  
   resetplt,/all
jump4:
; Compute mean direction weighted by energy
   mean_vec = [total(en*sin(phi)), total(en*cos(phi))]
   mean_dir = atan(mean_vec(0), mean_vec(1))*!radeg
   if mean_dir LT 0 THEN mean_dir = mean_dir + 360.	;my unders deg,from east
   return
END                             ; end angular_spec
;************************************************************
;PRO directions,up,vp,tp,rho,phi,ut,vt,uq,ui
PRO directions,up,vp,tp,rho,phi,ut,vt,uq,vq,ui,vi	;newly modified
; Procedure to compute mean direction of wave propagation
; from velocity and normalized temperature perturbations.
; Use covariances <u'T'+90> and <v'T'+90> to determine
; azimuth from North. Where T'+90 is the value of
; temperature after shifting phase by +90 deg via Hilbert
; transform.

   a = size(up)
   phi = fltarr(a(2))					;dir of wave propagation
   ut = phi & vt = phi
   uq = ut & ui = ut
   vq = ut						;newly added
   loop = a(2)-1					;number of profiles
   FOR k = 0, loop DO BEGIN
; Carry out Hilbert transform
      to = float(HILBERT(tp(*,k), -1))

; Low-pass filter the temperature to match smoothing
; inherent in the wind profiles (~300 m trangular weighted smoothing
; with Australian data)
      to = smooth(to, 5)
      ti = smooth(tp(*,k),5)    ; Untransformed temp
      x = total(to*up(*,k)*rho(*,k))/float(a(1))	;column average
      y = total(to*vp(*,k)*rho(*,k))/float(a(1))	;column average
      ut(k) = x
      vt(k) = y
      uq(k) = total(to*up(*,k))/n_elements(to)
      vq(k) = total(to*vp(*,k))/n_elements(to)		;newly added
      ui(k) = total(ti*up(*,k))/n_elements(ti)
; Phi is azimuth from north
      phi(k) = atan(x, y)
   ENDFOR
   return
END                             ; End Directions
;*********************************************************
PRO tspec_plot, ts, m,lamda
   a = size(ts)
   n = a(1)/2-10
   plot_oo,m(1:n),ts(1:n),xstyle = 9, $
    xmargin =[15,0],$
    xtitle='!17WAVENUMBER (cpm)',$
    ytitle='!17Normalized Temperature PSD (cpm!E-1!N)',/NODATA
   oplot,m(1:n),ts(1:n)
   axis,xaxis=1,xstyle = 1, $
    xrange=[lamda(1), $
            lamda(n)],xtitle='WAVELENGTH (km)'
   return
END
;*************************************************
PRO t_spec,t,tspec,err,m,lamda,navr,cpm
; Procedure to calulate vertical-wavenumber power spectrum for temperature
; m and lamda are the output vertical wavenumber and scale arrays
; averaged over navr points
; cpm is the sample spacing in metres
; tspec is the output array (real) of length npts/navr
; err is a standard deviation of the tspec at each wavenumber
; Generate window function
; Average over all available observations if necessary
   a=size(t)
   npts=a(1)
   IF (a(0) EQ 2) THEN n_obs = a(2) $
    ELSE n_obs = 1
   window = welch(npts)
   nred = npts/navr
   ts = fltarr(npts)
   tsq =fltarr(npts)
   nmax=npts-1
   nf = float(npts)/(total(window*window))

   delm = 1.0/(cpm*float(nred))
   m = findgen(nred)*delm
   if (navr ne 1) then m = m + delm/2.		; ? my question
   m(0) = delm/2.
   FOR i = 0,n_obs-1 DO BEGIN
      uf = complex((t(0:nmax,i)-total(t(0:nmax,i))/npts)*window ,0.0)
      ut = fft(uf,-1)
      spec = float(ut*conj(ut))
      ts = ts +spec
      tsq =spec*spec+tsq
   ENDFOR
   fnp =float(n_obs)
   norm = float(2*navr)*nf/(n_obs*delm) ; NB One-sided spectrum
   err = tsq/fnp - ts*ts/(fnp*fnp)
   err = fnp*sqrt(err)
   tspec = rebin(ts,nred)*norm
   tspec = tspec(0:npts/2)
   err = rebin(err,nred)*norm
   err = err(0:npts/2)
   lamda = 1.0/(m*1000)
return
END ; t_spec
;*********************************************************************
PRO rotary_spec,u,v,acp,ccp,m,lamda,navr,cpm
;Procedure to calulate vertical-wavenumber power spectrum for wind components
;m and lamda are the output vertical wavenumber and wavelength arrays
;averaged over navr points
;cpm is the sample spacing in metres
;uspec and vspec are the output arrays (real) of length npts/navr

;Generate window function
   a=size(u)
   npts=a(1)
   n2 = npts/2
   IF (a(0) EQ 2) THEN n_data = a(2) ELSE n_DATA =1
   window = welch(npts)
   nred = npts/navr
   nmax=npts-1
   nf = float(npts)/(total(window*window))
   delm = 1.0/(cpm*float(nred))
   m = findgen(nred)*delm
   IF (navr NE 1) THEN m = m + delm/2.
   m(0) = delm/2.
   norm = float(navr)*nf/delm
   acp = fltarr(n2) & ccp = acp

   FOR i = 0,n_data-1 DO BEGIN
      uf = (u(0:nmax,i)-total(u(0:nmax,i))/npts)*window
      vf = (v(0:nmax,i)-total(v(0:nmax,i))/npts)*window
      ut = fft(complex(uf,vf), -1)
      us = float(ut*conj(ut))
      IF (navr NE 1) THEN  $
       us = rebin(us,nred)*norm ELSE us = us*norm
;Anticlockwise component is positive half space
;clockwise component is negative half space
      acp = us(1:n2)+acp
      ccp = rotate(us(n2:nmax),2)+ccp
   ENDFOR
   acp = acp/float(n_data)
   ccp = ccp/float(n_data)
   lamda  = 1.0/(m*1000.0)
   return
END                             ; End rotary_spec
;*********************************************************************
;PRO constants,lat,n_bar,f,bo,co,wbar1,wbar2
PRO constants,lat,n_bar,f,bo,co,wbar1,wbar2,wbar1p,deltaw	;I add this
; Procedure to compute spectral constants
; based on Fritts and VanZandt, JAS, 1993.
; lat = station latitude in degrees
; n_bar = mean brunt (stability) frequency

   latt = abs(lat)
   inertial_period = 12./sin(latt*!dtor)
   f = 2.0*!pi/(3600.0*inertial_period)

   f_hat = f/n_bar
   p = 5./3.                    ; spectral slope of frequncy spectrum
   q = (p-1.)
   r = (2.-p)
   s = (p+1.)
   Bo = q*f^q/(1.-f_hat^q)
   Co = (f^(-q))*((1.-f_hat^q)/q - (1.-f_hat^s)/s)/(1.-f_hat^2)
; Compute mean frequency for canonical spectrum
; For p= 5/3
   p = 5./3.
   wbar1 = (q/r)*n_bar*(f_hat^q)*(1-f_hat^r)		;eq. (16) of VAE97
   wbar1p = (q/r)*n_bar*(f_hat^q)*(1-f_hat^(-r))	;I add this
   deltaw=1-f*f/wbar1p/wbar1p				;I add this
   wbar2 = -f*alog(f_hat)				;eq. (17) of VAE97
   return
END

;*********************************************************
PRO fit, ts,sigma,m,a,eo,ms,t
; Procedure to fit m-spectrum for normalized temperature spectrum (ts)
; to modified Desaubies spectrum in vertical wavenumber (m)

   dummystring = ' '
   npt = n_elements(m)          ; Make spectrum one-sided
   n = npt/2.-1
   s = ts(1:n)
   x = m(1:n)
   ss = sigma(1:n)
; Make initial guess of fitting constants
; Get mean-square normalized temperature
   dm = m(4)-m(3)
   ts_bar = total(s)*dm
; Next mean wavenumber
   m_bar = total(s*x)*dm/ts_bar
; First guess at m_star
   ms = x(1)
   fo = ts_bar/ms
   t = 4.0 
   a = [fo,ms,t ]
   fa = {X:x, Y:s, ERR:ss}
   functc,x,a,f        ; f is starting model 
   deviates = (s-f) 
; fitting a model to the avaerage spectrum using fitting algorithm
; embodied in mpfit.pro
   p_fit = mpfit('mp_funct',a,functargs= fa,/QUIET)
   functc,x,p_fit,spec_fit      ; evaluate fitted function
;loadct,40
;   Oplot,x,spec_fit,line = 1,color=240    ; plot it over data
;loadct,0
   a = p_fit
   return
END 
;created on 5/14/2002

PRO functc2,x,a,f,pder
;
; I add this
;
; Return a function value in f.
; based on lmfunct
; Function for fitting modified Desaubies spectrum
; 23 June 1998 RAV
; F = Ao mu/(1+mu^A2)
; where mu = x/A1

;
;a(0)=F_0
;a(1)=m_*
;a(2)=t+1
;
common fit6,t6_0
   mu = x/a(1)				;m/m*
;  d = mu^a(2)				;(m/m*)^(t+1)
   d = mu^(t6_0+1.0)
   b = 1./(1+d)
   c = mu*b
   f = a(0)*c				;f0 (m/m*) / (1+(m/m*)^(t+1))
   IF (n_params() GE 4) THEN $
    pder = [ [f/a(0)],  $
             [-f*b*(1+(1-(t6_0+1.0))*d)/a(1)] ]
   return
END 

FUNCTION mp_funct2, p,X=x,Y=y,ERR=sigma
;
; I add this
;
; this function provided to be called by the mpfit.pro
; routine . Devised Aug 12 1998 AC Beresford , Atmos Phys Univ of Adelaide
   nofx = n_elements(x)
   model = fltarr(nofx)
   deriv =fltarr(nofx)
;  functc, x,p,model
   functc2, x,p,model
   deriv = (y-model)/sigma
   return, deriv
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION Welch2,n
;similar to welch but with minor difference (Num Recipes in C, p 547)
;Procedure for calculating Welch Window function (Num Recipes, p 422)
j = findgen(n)
n2=n/2.
RETURN, 1. - ((j - n2)/n2)^2
END ; End Welch

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++
;--> nht (number of levels for a hrrd profile)
function get_eclipse_nht,fname

nht=0
blank = ' '
header = ' '

close,1
openr,1,fname
on_ioerror, break               ;IF I/O ERROR GO TO BREAK

for j = 1, 20 do readf,1,blank

;NOW READ IN DATA, SAVING DATA WHEN END OF FILE IS ENCOUNTERED
i = 0
repeat begin
    readf,1,header
    i=i+1
endrep until eof(1)
nht=i

break:   close,1
return,nht
end
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
pro gw_analysis2,up $			;zonal wind perturbation (m/s)
               ,vp $			;meridional wind perturbation (m/s)
	       ,arp $			;ascending rate perturbation (m/s)
               ,tp $			;temperature perturbation (K)
               ,um $			;mean zonal wind (m/s)
               ,vm $			;mean meridional wind (m/s)
	       ,arm $			;mean ascending rate (m/s)
               ,tb $			;mean temperature (K)
               ,pr $			;pressure (hpa)
               ,lat $			;latitude (deg)
               ,lon $			;longitude (deg)
               ,zkm $			;vertical grid (km)
               ,dz $			;vertical grid resolution (km)
               ,ipw $			;0 no prewhitening; 1 prewhitening
               ,iw $			;0 old method; 1 new method
               ,id $			;0 no linear detrend; 1 linear detrend
               ,ke,kez,kev,kevert,pe,f_w,m_bar,up_frac,mean_dir,uw,vw $ ;output
               ,n_bar=n_bar,plot=plot,ttl=ttl

common fit6,t6_0
nnn=n_elements(zkm)
;dz=zkm(1)-zkm(0)
dzm=dz*1000.
zm=zkm*1000.
z = zkm
p = pr*100.
p_r = p/287.05

tph = tp/tb
tph2=tph				;I add this
dt = fltarr(nnn)
dt = deriv(tb)/dz 			; Temperature gradient K/km
n_sq = (dt+9.8)*9.8/(tb*1000.)
n_bar = sqrt(avg1d(n_sq))

;print,'ok1'
;constants,lat,n_bar,f,bo,co,w1,w2
constants,lat,n_bar,f,bo,co,w1,w2,w1p,deltaw	;I add this
;print,'ok2'

; mean square value of normalized temperature perturbation
tps_bar =avg1d(tph^2)
;print,'ok3'

;linear detrend tph
if id eq 1 then begin
    tphp=tph
    detrend, tphp, z, 1
    tph2=tphp
endif

;prewhitening
if ipw eq 1 then begin
    tph2=tph(1:nnn-1)-tph(0:nnn-2)
endif

; Compute average spectrum
;print,'ok3a'
if iw eq 0 then t_spec_new,tph2,tsh,terr,m,lamda,1,dzm
;print,'ok3b'
if iw eq 1 then t_spec2,tph2,tsh,terr,m,lamda,1,dzm
;print,'ok3c'
;stop,'666'

; postdarkening
if ipw eq 1 then postdark, tsh, m

; Modify spectrum for response time
;tsh = tsh*(1+(40*m*2*!pi)^2)
; modify error estimates as well
; added 18 aug 1998 AC Beresford atmos phys Adelaide
;terr = terr*(1+(40*m*2*!pi)^2)

; Plot average spectrum
;if keyword_set(plot) then tspec_plot,tsh,m,lamda
if keyword_set(plot) then tspec_plot_new,tsh,m,lamda,ttl=ttl

 goto,jump_fit
; Fit average spectrum
; method one
;print,'ok3d1'
;save,file='test_fit.sav',tsh,terr,m,n_bar,tps_bar,bo,co
 fit,tsh,terr,m,A
;stop,'888'
;print,'ok3d2'
ms = a(1)
t = a(2)-1.
;print,'ok3e'

; method two
fit2,tsh,m,t2
; method three

fit3,tsh,terr,m,t2,A
ms3=a(1)
t3=a(2)-1.

; method four
fit4,tsh,m,t4

; method five
fit5,tsh,m,t5

; method six
t6_0=t5
fit6, tsh,terr,m,a
ms6 = a(1)

; method seven
t6_0=t
fit6, tsh,terr,m,a
ms7 = a(1)

; method eight
t6_0=t2
fit6, tsh,terr,m,a
ms8 = a(1)

; method nine
t6_0=t4
fit6, tsh,terr,m,a
ms9 = a(1)

; method ten
fit10,tsh,m,ms10,t10

; method eleven
t6_0=t5
fit11,tsh,m,ms11
jump_fit:

;print,'ok4'

; Compute wind perturbations.
; Average over central part of spectrum to avoid end effects
n1 = 5
n2 = nnn-6

upg = up(n1:n2)
vpg = vp(n1:n2)
arpg = arp(n1:n2)
umg = um(n1:n2)
vmg = vm(n1:n2)
armg = arm(n1:n2)
tpg = tp(n1:n2)/tb(n1:n2)
ke = 0.5*avg1d(upg^2+vpg^2)
pe = 0.5*avg1d((9.8*tpg/n_bar)^2)
en = ke+pe
;temp = moment(en,sdev = se)
et = en
tbg = tb(n1:n2)
kez = 0.5*avg1d(upg^2)
kev = 0.5*avg1d(vpg^2)
kevert = 0.5 * avg1d(arpg^2)
eo = tps_bar/(bo*co)*(9.8/n_bar)^2
;stop,'222'

; Compute density profile
p_rg = p_r(n1:n2)
rho=p_rg/tbg

;print,'ok5'

; Compute fraction of upgoing energy.
; In southern hemisphere upgoing waves have 
; anticlockwise polarisation (acp)
; In northern hemisphere use clockwise component (ccp).
rotary_spec,upg,vpg,acp,ccp,m,lamda,1,dz*1000.
IF (lat LT 0) THEN  $
    up_frac = total(acp)/total(acp+ccp) $
ELSE  $
    up_frac = total(ccp)/total(acp+ccp)

;print,'ok6'

; Compute mean direction of propagation
;directions,upg,vpg,tpg,rho,phi,ut,vt,uq,ui
directions_new,upg,vpg,tpg,rho,phi,ut,vt,uq,vq,ui,vi	;newly modified

;plot the angular distribution and calculate mean_dir
dir = phi*!radeg
IF dir lt 0 THEN dir=dir+360.
mean_dir=dir

;print,'ok7'
; Compute density-weighted momentum fluxes 
; following Vincent et al (1996)
fac = 9.8/n_bar^2
uw = fac*w1*ut			;energy weighted flux
vw = fac*w1*vt			;energy weighted flux
uw2 = fac*w1*uq			;newly added
vw2 = fac*w1*vq			;newly added
uwp = -fac*w1p*ut*deltaw	;newly added
vwp = -fac*w1p*vt*deltaw	;newly added
uw2p = -fac*w1p*uq*deltaw	;newly added
vw2p = -fac*w1p*vq*deltaw	;newly added

;print,'ok8'
; Compute mean frequency
stokes_new,upg,vpg,umg,vmg,phi,n_bar,f_w,df,zm(n1:n2)
;print,'ok9'

if iw eq 0 then phase_speed_new,umg,vmg,phi,upg,vpg, $
    f,n_bar,f_w,m_bar,k_bar,omega,c_z,c_i,c_x,c_y,df, $
    c_ix,c_iy,c_gx,c_gy,c_xcomp,c_ycomp, u_mean, u_mm, v_mm, theta,dzm
if iw eq 1 then phase_speed2,umg,vmg,phi,upg,vpg, $
    f,n_bar,f_w,m_bar,k_bar,omega,c_z,c_i,c_x,c_y,df, $
    c_ix,c_iy,c_gx,c_gy,c_xcomp,c_ycomp, u_mean, u_mm, v_mm, theta

;print,'ok10'
;save,filename=outflname, $
;    w1,w2,f,ke,pe,phi,ut,vt,f_w,df,n_bar,m_bar,k_bar, $
;    omega,c_z,c_i,c_x,c_y,c_ix,c_iy,c_gx,c_gy,kez,kev, $
;    sez,en,se,c_xcomp,c_ycomp,month,et,lat,lon,eo, $
;    up_frac,uw,vw,mean_dir,uq,vq,ui,acp,ccp,type,ms,t,u_mean, $
;    theta,u_mm,v_mm, $
;    uw2,vw2,notes,t2,ms3,t3,t4,t5,ms6,ms7,ms8,ms9,ms10,t10,ms11, $
;    uwp,vwp,uw2p,vw2p
;save,filename=outflname,ke,kez,kev,pe,f_w,m_bar,up_frac,mean_dir,uw,vw,notes,k_bar,lat
;print,'ok11'

return
end

;CALCULATE THE COROLIOS PARAMETER FROM LATITUDE
;APPLICABLE TO BOTH SCALAR AND ARRAY OF LATITUDES
function get_coriolis, latitude ;LATITUDE, IN DEG
omega = 7.2919996d-05
f = 2.0d*omega*sin(latitude*!dtor)
return, f
end

;( INTRINSIC FREQUENCY, VERTICAL WAVENUMBER) ==> HORIZONTAL WAVENUMBER

function dispersion_m_w_to_k,w,m,lat $
                            ,n=n,hrou=hrou,help=help
;----------------------------------------------
if keyword_set(help) then begin
    print,''
    print,'Purpose: '
    print,'         ( intrinsic frequency + vertical wavenumber)'
    print,'                 --> horizontal wavenumber'
    print,'Usage: '
    print,'         k=dispersion_m_w_to_k(w,m,f,n=n,alpha=alpha)'
    print,'Steps:'
    print,'         '
    print,'Input:'
    print,'         w -- intrinsic frequency [SI]'
    print,'         m -- vertical wavenumber [SI], 2*!PI / Lz'
    print,'         lat -- latitude [deg]'
    print,'Keywords:'
    print,'         help -- print this information'
    print,'         n -- Brunt-Vasalla freq [SI]'
    print,'         hrou -- density scale height [km]'
    print,'Output:'
    print,'         k -- horizontal wavenumber [SI], 2*!PI / Lh'
    print,'Examples:'
    print,'         '
    print,'History:'
    print,'         08/11/2006 created by Jie Gong'
    print,'Note:'
    print,'         '
    print,''
    stop
endif
;----------------------------------------------
if ~ keyword_set(hrou) then hrou=7.
if ~ keyword_set(n) then n=0.02                 ;typical stratospheric value
alpha=0.5/(hrou*1.e3)
f=get_coriolis(lat)
k=sqrt((w*w-f*f)/(n*n-w*w)*(m*m+alpha*alpha))
;----------------------------------------------
return,k
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;c

function gw_analysis_struct,input,plot=plot,ttl=ttl,help=help
;----------------------------------------------
if keyword_set(help) then begin
    print,''
    print,'Purpose: '
    print,'         do the GW analysis following VAE97'
    print,'Usage: '
    print,'         gw=gw_analysis_struct(input)'
    print,'Steps:'
    print,'         '
    print,'Input:'
    print,'         input -- the input structure'
    print,'             .up : zonal wind perturbation [m/s]'
    print,'             .vp : meridional wind perturbation [m/s]'
    print,'		.arp: ascending rate perturbation [m/s]'
    print,'             .tp : temperature perturbation [K]'
    print,'             .um : mean zonal wind [m/s]'
    print,'             .vm : mean meridional wind [m/s]'
    print,'		.arm: mean ascending rate [m/s]'
    print,'             .tm : mean temperature [K]'
    print,'             .pm : mean pressure [hpa]'
    print,'             .lat : latitude location [deg]'
    print,'             .lon : longitude location [deg]'
    print,'             .z : vertical grid [km]'
    print,'             .dz : vertical grid resolution [km]'
    print,'             .ipw : 0 no prewhitening; 1 prewhitening'
    print,'             .iw : 0 old method; 1 new method'
    print,'             .id : 0 no linear detrend; 1 linear detrend'
    print,'Keywords:'
    print,'         help -- print this information'
    print,'         plot -- make preliminary plot'
    print,'         ttl -- title of the preliminary plot'
    print,'Output:'
    print,'         gw -- all the GW info. derived from this analysis'
    print,'               see my dissertation or VAE97 for more details'
    print,'             .ke : kinietic energy density [SI]'
    print,'             .kez : zonal kinietic energy density [SI]'
    print,'             .kev : meridional kinietic energy density [SI]'
    print,'             .pe : potential energy density [SI]'
    print,'             .f_w : intrinsic frequency/f'
    print,'             .m_bar : reciprocal of vertical wavelength [SI]'
    print,'             .up_frac: fraction of upward propagating waves'
    print,'             .mean_dir : mean propagation direction [deg]'
    print,'                         azimuth from north (?)'
    print,'             .uw : zonal flux [?]'
    print,'             .vw : meridional flux [?]'
    print,'             .n_bar : buoyancy frequency [SI]'
    print,'             .lz: vertical wavelength [km]'
    print,'             .lh: horizontal wavelength [km]'
    print,'             .m: vertical wavenumber, m=2*Pi*m_bar [SI]'
    print,'             .f: Coriolis parameter [SI]'
    print,'             .lat: latitude [deg]'
    print,'             .omega: intrinsic frequency [SI]'
    print,'             .k: horizontal wavenumber [SI]'
    print,'             .kx: zonal wavenumber [SI]'
    print,'             .ky: meridional wavenumber [SI]'
    print,'             .input: input to calculate the gw para.'
    print,'		.kevert: vertical kinetic energy density [SI]'
    print,'Note:'
    print,'         virtually the same as gw_analysis2, except to make the'
    print,'         code more reader friendly and easier to use'
    print,''
    stop
endif
;----------------------------------------------
gw_analysis2,input.up,input.vp,input.arp $
	    ,input.tp,input.um,input.vm,input.arm $
            ,input.tm,input.pm,input.lat[0],input.lon[0] $
            ,input.z,input.dz,input.ipw,input.iw,input.id $
            ,ke,kez,kev,kevert,pe,f_w,m_bar,up_frac,mean_dir,uw,vw $
            ,n_bar=n_bar,plot=plot,ttl=ttl
;----------------------------------------------
lz=1.e-3/m_bar					;[km]
m=2.*!pi*m_bar					;2*!PI / Lz instead of 1/Lz
lat=input.lat[0]					;[deg]
f=get_coriolis(lat)
omega=f_w*f					;intrinsic frequency [SI]

k=dispersion_m_w_to_k(omega,m,lat,n=n_bar,hrou=7.)
lh=2.e-3*!pi/k					;horizontal wavelength [km]
phi2kl,mean_dir,lh,kx,ky
gw={ke:ke,kez:kez,kev:kev,kevert:kevert,pe:pe,f_w:f_w,m_bar:m_bar,$
    up_frac:up_frac, $
    mean_dir:mean_dir,uw:uw,vw:vw,n_bar:n_bar, $
    lz:lz,lh:lh,m:m,f:f,lat:lat,omega:omega,k:k,kx:kx,ky:ky,input:input,$
    outputflag:0}
;----------------------------------------------
return,gw
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; PURPOSE : cal. the perturbations using a selected method for one sounding
;
; INPUT :
;
; OUTPUT :
;
; ALGORITHM :
;
; NOTES : created on 1/26/2004
;
;         modified on 4/6/04 to make the code more efficient
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro perturb_hrrd $

               ;input
               ,u $			;u interpolated to regular grid
               ,flag $			;0: enough valid data; 1: not enough
               ,npoly=npoly $		;order of polynomial fit
               ,flow=flow $		;lower limit of the high-pass filter
               ,fhigh=fhigh $		;upper limit of the low-pass filter

               ;output
               ,up $			;perturbation, same unit as the input
               ,ub			;background, same unit as the input

if flag eq 1 then return

;polynomial fit
if keyword_set(npoly) then begin
    co=poly_fit(double(findgen(n_elements(u))),double(u),npoly,ub)
    up=u-ub
    return
endif
;high-pass filter --> perturbation
if keyword_set(flow) then begin
    up=digital_smooth(u,flow,1)
    ub=u-up
    return
endif
;low-pass filter --> background
if keyword_set(fhigh) then begin
    ub=digital_smooth(u,0,fhigh)
    up=u-ub
    return
endif

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; PURPOSE : obtain the raw data info. (data, missing, code) for one radiosonde
;
; INPUT : name of radiosonde profile (fname)
;
; OUTPUT : data,missing,code
;
; ALGORIGHM :
;
; NOTES : more complete than read_whole_profile.pro
;
; missing data flag: 	P	mb	9999.0
;			T	c	999.0
;			U	m/s	9999.0
;			V	m/s	9999.0
;			AR	m/s	999.0
;			Z	m	99999.0
;          see the more complete in missing
;
;HISTORY: change read_hrrd to read_eclipse to handle eclipse balloon sounding
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;c

pro read_eclipse $

    ;INPUT
    ,fname $			;name of the Sounding
    ,fast=fast $		;if set, faster but larger strcuture
				;also there may be zero paddings
				;so it's better not to use this option

    ;OUTPUT
    ,data $			;data
    ,nrec $			;record length
    ,missing $			;missing data info.
    ,code			;quality control info.

;--> nht
if ~ keyword_set(fast) then begin
    nht=get_eclipse_nht(fname)-10 ;with newly processed data, the last 10 lines dedicated to others
endif else begin
    nht=8000			;MAX POSSIBLE NUMBER OF OBS. FOR ONE PROFILE
endelse

data = {tm:fltarr(nht), $	;Time			s	9999.0
        p:fltarr(nht), $	;Pressure		hPa	9999.0
        t:fltarr(nht), $	;Dry-bulb temp		C	999.0
        rh:fltarr(nht), $	;Relative Humidity	%	999.0
        ws:fltarr(nht), $	;Wind Speed		m/s	999.0
        dir:fltarr(nht), $	;Wind Direction		deg	999.0
        lon:fltarr(nht), $	;Longitude		deg	999.0
        lat:fltarr(nht), $	;Latitude		deg	999.0
        z:fltarr(nht), $	;Altitude		m	99999.0
	gph:fltarr(nht),$	;Geopotential Hgt	m	99999.0
	mri:fltarr(nht),$	;??			C	999.0	
	ri:fltarr(nht),$	;??			C	999.0
        tdew:fltarr(nht), $	;Dew Point		C	999.0
        tvir:fltarr(nht), $	;Virt. temp		C	999.0
	rs:fltarr(nht),$	;Rising speed		m/s	999.0
        ele:fltarr(nht), $	;Elevation Angle	deg	999.0
        azi:fltarr(nht), $	;Azimuth		deg	999.0
	range:fltarr(nht),$	;			m	999.0
	D:fltarr(nht),$		;			kg/m3 	999.0
	u:fltarr(nht),$		;u wind,derived 	m/s	999.0
	v:fltarr(nht)$		;u wind,derived 	m/s	999.0
	}
missing = {tm:9999.0, $
           p:9999.0, $
           t:999.0, $
           rh:999.0, $
           ws:999.0, $
           dir:999.0, $
           lon:999.0, $
           lat:999.0, $
           z:99999.0, $
           gph:99999.0, $
           mri:999.0, $
           ri:999.0, $
           tdew:999.0, $
           tvir:999.0, $
           rs:999.0, $
           ele:999.0, $
           azi:999.0, $
           range:999.0, $
           D:999.0, $
	   u:999.0, $
	   v:999.0 $
           }
;code='99.0    Unchecked;'+ $
;     '1.0     Checked (Good);'+ $
;     '2.0     Checked (Maybe);'+ $
;     '3.0     Checked (Bad);'+ $
;     '4.0     Checked (Estimated);'+ $
;     '9.0     Checked (Missing)'
;

blank = ' '
header = ' '

close,1
openr,1,fname
on_ioerror, break		;IF I/O ERROR GO TO BREAK

for j = 1, 20 do readf,1,blank	

;NOW READ IN DATA, SAVING DATA WHEN END OF FILE IS ENCOUNTERED
i = 0
iii = 0
repeat begin
    readf,1,header
    iii=strsplit(header,count=c)
    if c NE 19 then goto,stop_reading
    data.tm(i)=float(strmid(header,iii(0),iii(1)-iii(0)))
    data.p(i)=float(strmid(header,iii(1),iii(2)-iii(1)))
    data.t(i)=float(strmid(header,iii(2),iii(3)-iii(2)))
    data.rh(i)=float(strmid(header,iii(3),iii(4)-iii(3)))
    data.ws(i)=float(strmid(header,iii(4),iii(5)-iii(4)))
    data.dir(i)=float(strmid(header,iii(5),iii(6)-iii(5)))
    data.lon(i)=float(strmid(header,iii(6),iii(7)-iii(6)))
    data.lat(i)=float(strmid(header,iii(7),iii(8)-iii(7)))
    data.z(i)=float(strmid(header,iii(8),iii(9)-iii(8)))
    data.gph(i)=float(strmid(header,iii(9),iii(10)-iii(9)))
    data.mri(i)=float(strmid(header,iii(10),iii(11)-iii(10)))
    data.ri(i)=float(strmid(header,iii(11),iii(12)-iii(11)))
    data.tdew(i)=float(strmid(header,iii(12),iii(13)-iii(12)))
    data.tvir(i)=float(strmid(header,iii(13),iii(14)-iii(13)))
    data.rs(i)=float(strmid(header,iii(14),iii(15)-iii(14)))
    data.ele(i)=float(strmid(header,iii(15),iii(16)-iii(15)))
    data.azi(i)=float(strmid(header,iii(16),iii(17)-iii(16)))
    data.range(i)=float(strmid(header,iii(17),iii(18)-iii(17)))
    data.D(i)=float(strmid(header,iii(18),strlen(header)-iii(18)))
    i=i+1

endrep until eof(1)

break: close,1
stop_reading:
nrec=i

data.u(0:nrec-1)=data.ws(0:nrec-1)*sin(data.dir(0:nrec-1)/180.*!pi)
data.v(0:nrec-1)=data.ws(0:nrec-1)*cos(data.dir(0:nrec-1)/180.*!pi)

return
end
;CREATE AN MONOTONIC INCREASING ARRAY OF increment dy, WITH MIN Y0, MAX Y1
function make_array3, y0, y1, dy $
                    ,loose=loose $	;not strict about integral multiples
                    ,count=count
on_error,2
if y0 ge y1 then begin
    if y0 gt y1 then message,'Error: y0 should be no bigger than y1' else y=y0
endif else begin
    n=long((y1-y0)/float(dy)+1.1)
    y = fltarr(n)
    y = y0+findgen(n)*dy
    if (abs(y(n-1)-y1) gt 0.1*dy) and $
       (not keyword_set(loose)) then message,'Error: dy not intergral increment'
endelse
count=n_elements(y)
return, y
end
pro check_hrrd_gap $

                  ;input
                  ,zraw $	;raw z data
                  ,uraw $	;raw u data (input/output)
                  ,m_z $	;missing flag for z
                  ,m_u $	;missing flag for u
                  ,z0 $		;bottom z (m)
                  ,z1 $		;top z (m)
                  ,z0p $	;extended z0 (m)
                  ,z1p $	;extended z1 (m)
                  ,note $	;descriptive note on the file being analyzed
                  ,madz=madz $	;if set, the maximum allowed valud of gap (m)
				;default=1 km
                  ,bndy=bndy $	;if set, the maximum gap in the boundary (m)
				;default=0 km

                  ;output
                  ,flag $	;0: enough data; 1: not enough data
                  ,z_out $ 	;selected z
                  ,u_out 	;selected u

;(0) predefine some variables
flag=0	;0: enough data, ub,up successfully retrieved ;1: not enough data
if not keyword_set(madz) then madz=1000.
if not keyword_set(bndy) then bndy=0.
z=zraw
u=uraw

;(1) pick up z,u that correspond to where z and u are not missing
ava=where(z ne m_z and u ne m_u,nava)
if nava eq 0 then begin
    print,note+' No data at all'
    flag=1
    return
endif
z=z(ava)
u=u(ava)
;(2) pick up z,u where z is within the extended range [z0p,z1p]
ava=where(z ge z0p and z le z1p,nava)
if nava eq 0 then begin
    print,note+' No data within the extended altitude range'
    flag=1
    return
endif
z=z(ava)
u=u(ava)
;(3) check the boundaries, if ztop .le. z1-1km or zbot .ge. z0+1km, return
if max(z) le z1-bndy or min(z) ge z0+bndy then begin
    print,note+' Not enough boundary coverage'
    flag=1
    return
endif
;(4) sort z so that z is in ascending order, excluding same level data
ava=uniq(z,sort(z))
nz=n_elements(ava)
nz0=n_elements(z)
;if nz ne nz0 then begin
;    print,note+' has '+string2(nz0-nz)+' same-z data removed'
;endif
z=z(ava)
u=u(ava)
;(5) check the data gap inside, if larger than 1 km, return
ddzz=z(1:nz-1)-z(0:nz-2)
if max(ddzz) ge madz then begin
    print,note+' Too large data gap inside'
    flag=1
    return
endif

u_out=u
z_out=z
return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; PURPOSE : interpolate the hrrd data into regular grid
;
; INPUT :
;
; OUTPUT :
;
; ALGORITHM :
;
; NOTES : created on 1/26/2004
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro regular_hrrd $

    ;input
    ,zbig $			;common extended z-grid
    ,z $			;raw z in selected extended altitude range
    ,u $			;raw u in selected extended altitude range
    ,flag $			;enough data available/not available

    ;output
    ,ugrid $			;u interpolated to regular grid
    ,zgrid			;corresponding regular z-grid

nodata=-999.

;do not have enough valid data
if flag eq 1 then begin
    ugrid=nodata
    zgrid=nodata
    return
endif

;--> zgrid, the extended z-grid (m) specific to the given variable
i0=(where(zbig ge min(z)))(0)
i1=(where(zbig ge max(z)))(0)-1
;zgrid=zbig(projection(zbig,z,/irregular))
;zgrid=zbig(projection(zbig,z))
zgrid=zbig(i0:i1)

; interpolate to the extended z-grid
ugrid=spline(z,u,zgrid)

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; PURPOSE : map [u,v,t,p][p,b] in zgrid[[u,v,t,p] to the desired zkm
;
; INPUT :
;
; OUTPUT :
;
; ALGORITHM :
;
; NOTES : create on 1/26/2004
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro map_hrrd $

            ;input
            ,flag $	;0: enough valid data; 1: not enough
            ,zgrid $	;extended z-grid (m) specific to the given variable
            ,zm $	;common z-grid (desired)

            ;input/output
            ,up $	;perturbation
            ,ub		;background

nz=n_elements(zm)
;not enough data
if flag eq 1 then begin
    up=replicate(-999.,nz)
    ub=replicate(-999.,nz)
    return
endif else begin
    i0=(where(zgrid ge min(zm)))(0)
    i1=(where(zgrid ge max(zm)))(0)-1
    up=up(i0:i1)
    ub=ub(i0:i1)
endelse

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; PURPOSE : test and plot the perturbation/background for one hrrd sounding
;
; INPUT :
;
; OUTPUT :
;
; ALGORITHM :
;
; NOTES : created on 1/26/2004
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro test_hrrd_perturb,z,u,z2,u2,zgrid,ugrid,zkm,ub,up,flag,m_z,m_u,z0p,z1p $
                     ,fit=fit $			;if set, test the interpolation
                     ,perturb=perturb $		;if set, test the perturbation
                     ,id=id 			;if set, plot the id

if not keyword_set(id) then id=''
ava=where(z ne m_z and u ne m_u,nava)
if nava gt 0 then begin
    if keyword_set(fit) then begin
        loadct,0
        plot,u(ava),z(ava)/1000.,yran=[z0p/1000.,z1p/1000.],xst=1,yst=1,thick=4
        if flag eq 0 then begin
            oplot,u2,z2/1000.,linestyle=1,thick=6
            loadct,40
            oplot,u2,z2/1000.,thick=1,color=240
            loadct,34
            oplot,ugrid,zgrid/1000.,linestyle=2,thick=2,color=74
            loadct,39
            oplot,ub,zkm,linestyle=2,thick=2,color=120
        endif
    endif
    if keyword_set(perturb) and flag eq 0 then begin
        loadct,0
        plot,up,zkm,xst=1,yst=1,thick=4,title=id
    endif
endif

return
end

function digital_smooth,raw_data0,flow,fhigh
;  a low-pass filter is applied with flow=0,fhigh=0.1,nterm=20

   n=n_elements(raw_data0)
   raw_data=[replicate(raw_data0(0),10),raw_data0,replicate(raw_data0(n-1),10)]
   sz = size(raw_data)
   dim = sz(1)
   nterm = 2*n_elements(raw_data0)
   data_in = fltarr(dim+2*nterm)
   data_in(0:nterm-1) = raw_data(0)
   data_in(nterm+dim:dim+2*nterm-1) = raw_data(n_elements(raw_data)-1)
   data_in(nterm:nterm+dim-1) = raw_data
   a = 50
   coeff = digital_filter(double(flow),double(fhigh),double(a),nterm)
;  coeff = digital_filter((flow),(fhigh),(a),nterm)
   data_out = convol(double(data_in),double(coeff))
;  data_out = convol((data_in),(coeff))
   ;as well as the zero padding elements, the end 5 points are discarded
   smooth_data = float(data_out(nterm+10:nterm+dim-1-10))
;  smooth_data = float(data_out(nterm:nterm+dim-1))
   ;plot,raw_data(5:dim-6)
   ;plot,smooth_data
   ;stop
   return,smooth_data
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; PURPOSE : analyze one hrrd
;
; INPUT : fname                  ;name of single radiosonde sounding file
;         z0                     ;bottom z (m)
;         z1                     ;top z (m)
;         buffer                 ;buffer*dz = extended range, e.g., 34
;         dz                     ;vertical resolution (m)
;         lat                    ;latitude (deg)
;         outname		 ;file name of the output
;
; OUTPUT : create a file of name outname which contains the 
;
; ALGORITHM : [1] read raw data
;             [2] select the right raw data
;             [3] interpolate the raw data into regular grid
;             [4] get the perturbations and background fields
;             [5] do gw analysis
;
; EXAMPLES : gw_eclipse_complete,'../data/eclipse/T36_0230_121520_Artemis_Rerun_CLEAN.txt',$
;		2000.,8900.,3,30.,-39.2,gw=gw,npoly=2
;            3 is buffer;30.is dz,45.is lat.
;	     the major results are stored in the structure "gw"
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro gw_eclipse_complete $

               ;all are input
               ,fname $			;name of hrrd
               ,z0 $			;bottom z (m)
               ,z1 $			;top z (m)
               ,buffer $		;buffer*dz = extended range, e.g., 34
               ,dz $			;vertical resolution (m)
               ,lat $			;latitude (deg)
               ,cliptop=cliptop $	;if set, do not extend the top
               ,clipbot=clipbot $	;if set, do not extend the bottom
               ,draw=draw $		;if set, plot the perturbations, etc.
					;obsolete, use uplot etc
               ,uplot=uplot $		;plot the perturbations, etc. for u
               ,vplot=vplot $		;plot the perturbations, etc. for v
               ,tplot=tplot $		;plot the perturbations, etc. for t
               ,pplot=pplot $		;plot the perturbations, etc. for p
               ,nogw=nogw $		;if set, do not analyze gw
               ,npoly=npoly $		;order of polynomial fit
               ,flow=flow $ 		;lower limit of the high-pass filter
               ,fhigh=fhigh $ 		;upper limit of the low-pass filter
               ,id=id $ 		;id+analysis method of the sounding
               ,gw=gw $ 		;id+analysis method of the sounding
               ,notes=notes  		;notes of code

if ~ keyword_set(id) then id=''
if ~ keyword_set(notes) then notes=''

;[A] read in raw data and setting up some universal parameters
read_eclipse,fname,data,nrec,missing,code

z=data.z(0:nrec-1)
u=data.u(0:nrec-1)
v=data.v(0:nrec-1)
ar=data.rs(0:nrec-1)
p=data.p(0:nrec-1)
t=data.t(0:nrec-1)
lat=data.lat(0:nrec-1)
lon=data.lon(0:nrec-1)
m_z=missing.z
m_u=missing.u
m_v=missing.v
m_ar=missing.rs
m_p=missing.p
m_t=missing.t
m_lat=missing.lat
m_lon=missing.lon

if ~ keyword_set(clipbot) then z0p=z0-buffer*dz else z0p=z0-5*dz
if ~ keyword_set(cliptop) then z1p=z1+buffer*dz else z1p=z1+5*dz
zm=make_array3(z0,z1,dz)		;z-grid (m)
zkm=zm/1000.				;z-grid (km)
zbig=make_array3(z0p,z1p,dz)		;the common extended z-grid (m)

;[B] select the approximate data
;--> flag[u,v,t,p]	:		enough data available/not available
;    [u,v,t,p]2		:		raw data in the extended z-range
;    z2[u,v,t,p]	:		the corresponding z

inputflag=0
check_hrrd_gap,z,u,m_z,m_u,z0,z1,z0p,z1p,'U',flagu,z2u,u2
check_hrrd_gap,z,v,m_z,m_v,z0,z1,z0p,z1p,'V',flagv,z2v,v2
check_hrrd_gap,z,ar,m_z,m_ar,z0,z1,z0p,z1p,'AR',flagar,z2ar,ar2
check_hrrd_gap,z,t,m_z,m_t,z0,z1,z0p,z1p,'T',flagt,z2t,t2
check_hrrd_gap,z,p,m_z,m_p,z0,z1,z0p,z1p,'P',flagp,z2p,p2

if flagu EQ 1 or flagv EQ 1 or flagar EQ 1 or flagt EQ 1 $
	or flagp EQ 1 then begin
    inputflag=1
    goto,output_round 
endif
check_hrrd_gap,z,lat,m_z,m_lat,z0,z1,z0p,z1p,'Lat',flaglat,z2lat,lat2
check_hrrd_gap,z,lon,m_z,m_lon,z0,z1,z0p,z1p,'Lon',flaglon,z2lon,lon2

;[C] interpolate the data into regular grid
;--> [u,v,t,p]grid	:		raw data interpolated to regular grid
;    zgrid[u,v,t,p]	:		the corresponding z-grid

regular_hrrd,zbig,z2u,u2,flagu,ugrid,zgridu
regular_hrrd,zbig,z2v,v2,flagv,vgrid,zgridv
regular_hrrd,zbig,z2ar,ar2,flagar,argrid,zgridar
regular_hrrd,zbig,z2t,t2,flagt,tgrid,zgridt
regular_hrrd,zbig,z2p,p2,flagp,pgrid,zgridp
regular_hrrd,zbig,z2lat,lat2,flaglat,latgrid,zgridlat
regular_hrrd,zbig,z2lon,lon2,flaglon,longrid,zgridlon

;[D] cal. the perturbations at zgrid[u,v,t,p]
;--> [u,v,t,p][p,b]
perturb_hrrd,ugrid,flagu,npoly=npoly,flow=flow,fhigh=fhigh,up,ub
perturb_hrrd,vgrid,flagv,npoly=npoly,flow=flow,fhigh=fhigh,vp,vb
perturb_hrrd,argrid,flagar,npoly=npoly,flow=flow,fhigh=fhigh,arp,arb
perturb_hrrd,tgrid,flagt,npoly=npoly,flow=flow,fhigh=fhigh,tp,tb
perturb_hrrd,pgrid,flagp,npoly=npoly,flow=flow,fhigh=fhigh,pp,pb
perturb_hrrd,latgrid,flaglat,npoly=npoly,flow=flow,fhigh=fhigh,latp,latb
perturb_hrrd,longrid,flaglon,npoly=npoly,flow=flow,fhigh=fhigh,lonp,lonb


;[E] map [u,v,t,p][p,b] in zgrid[[u,v,t,p] to the desired zm/zkm
map_hrrd,flagu,zgridu,zm,up,ub
map_hrrd,flagv,zgridv,zm,vp,vb
map_hrrd,flagar,zgridar,zm,arp,arb
map_hrrd,flagt,zgridt,zm,tp,tb
map_hrrd,flagp,zgridp,zm,pp,pb
map_hrrd,flaglat,zgridlat,zm,latp,latb
map_hrrd,flaglon,zgridlon,zm,lonp,lonb

;[F] test and plot the perturbation/background
if keyword_set(draw) then $
    test_hrrd_perturb,z,u,z2u,u2,zgridu,ugrid,zkm,ub,up,flagu,m_z,m_u,z0p,z1p $
                     ,/perturb,id=id+' U'
if keyword_set(uplot) then $
    test_hrrd_perturb,z,u,z2u,u2,zgridu,ugrid,zkm,ub,up,flagu,m_z,m_u,z0p,z1p $
                     ,/perturb,id=id+' U'
if keyword_set(vplot) then $
    test_hrrd_perturb,z,v,z2v,v2,zgridv,vgrid,zkm,vb,vp,flagv,m_z,m_v,z0p,z1p $
                     ,/perturb,id=id+' V'
if keyword_set(vplot) then $
    test_hrrd_perturb,z,v,z2ar,ar2,zgridar,argrid,zkm,arb,arp,flagar,m_z,m_ar,$
    z0p,z1p,/perturb,id=id+' AR'
if keyword_set(wplot) then $
    test_hrrd_perturb,z,t,z2t,t2,zgridt,tgrid,zkm,tb,tp,flagt,m_z,m_t,z0p,z1p $
                     ,/perturb,id=id+' T'
if keyword_set(pplot) then $
    test_hrrd_perturb,z,p,z2p,p2,zgridp,pgrid,zkm,pb,pp,flagp,m_z,m_p,z0p,z1p $
                     ,/perturb,id=id+' P'

;[G] convert the unit of T
if flagt eq 0 then tb=tb+273.16		;C-->K

output_round:
;[H] do the gravity wave analysis
if ~ keyword_set(nogw) and inputflag EQ 0 then begin
    input={up:up,vp:vp,arp:arp,tp:tp,um:ub,vm:vb,arm:arb,$
	    tm:tb,pm:pb,lat:latb,lon:lonb,z:zm*1.e-3,dz:dz*1.e-3,ipw:0,iw:0,id:0}
    gw=gw_analysis_struct(input)
endif
if keyword_set(nogw) or inputflag NE 0 then begin
    gw={ke:-999.,outputflag:1}
endif

return
end
