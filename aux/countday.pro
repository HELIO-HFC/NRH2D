;+
;NAME:
;		countday
;
;PURPOSE:
;		Returns a string vector which contains the dates between 
;		the two input days. 
;
;CALLING SEQUENCE:
;		IDL>Result = countday(first_day,last_day,nday=nday,julian=julian,time=time)
;
;INPUTS:
;		first_day - Scalar containing the first day to return.
;					(The format can be YYYYMMDD or YYYY-MM-DD, where YYYY is the year,
;					 MM the month, and DD the day of the first day.)
;		last_day  - Scalar containing the last day to return.
;					(The format can be YYYYMMDD or YYYY-MM-DD, where YYYY is the year,
;					 MM the month, and DD the day of the first day.)
;		
;OPTIONAL INPUTS:
;		nday	- Specifies the number of day to return.
;				  If nday is given, then the value of the input parameter 
;				  last_day is ignored, and this latter is considered
;				  as an output parameter which returns the last date of the
;				  output vector, last_day = first_day + nday.  
;		time	- Specifies the time of the day. 
;				  (only useful for the output parameter julian.)
;
;KEYWORD PARAMETERS:
;		None.
;
;OUTPUTS:
;		days	- String vector containing the dates between first_day and last_day.
;				  (or between first_day and first_day + nday if nday is specified.)
;
;OPTIONAL OUTPUTS:
;		nday	- If nday is not specify as an optional input parameter,
;				  then nday returns the number of days between first_day and last_day.
;				  (i.e. the number of elements of the output dates.)
;		julian	- Returns the corresponding julian days.
;
;EXAMPLE:
;		IDL>Days = countday('20010101','20010104',nday=nday)
;		IDL>for i=0,nday-1 do print,Days[i]
;		20010101
;		20010102
;		20010103
;		20010104
;
;MODIFICATION HISTORY:
;	28-FEB-2006,X.Bonnin:	Written.
;	04-MAR-2007,X.Bonnin:	Added the optional input parameter time.
;	24-OCT-2007,X.Bonnin:	Renamed input parameter day to first_day.
;	26-APR-2010,X.Bonnin:	last_day was now an input parameter instead of
;							a optional input parameter.
;	20-DEC-2011,X.Bonnin:	Input format YYYY-MM-DD is also correct.
;-

FUNCTION countday,first_day,last_day, $
				  nday=nday, $
				  time=time,julian=julian

;Test input parameters
if (~keyword_set(first_day)) then begin
    message,/info,'Call is:'
    print,'Result = countday(first_day,last_day, $'
    print,'                  nday=nday, $'
    print,'                  time=time,julian=julian)'
    return,''
endif

;Split the input date
fday = strsplit(strtrim(first_day[0],2),'-',/EXTRACT)
fflag = n_elements(fday)
case fflag of
	1:begin
		yyyy = strmid(fday,0,4) & mm = strmid(fday,4,2) & dd = strmid(fday,6,2)
	end
	3:begin
		yyyy = fday(0) & mm = fday(1) & dd = fday(2)
	end
	else:begin
		message,/CONT,'Wrong input format!'
		return,''
	end
end

;Time 
if (~keyword_set(time)) then time = '12:00:00' else t = strtrim(time[0],2) 
hrs = strsplit(time,':',/EXTRACT)
hh = hrs[0] & mn = hrs[1] & ss = hrs[2]
hrs = 0B ;free memory

;Get the Julian day for fday
fjul = (julday(mm,dd,yyyy,hh,mm,ss))[0]


if (keyword_set(nday)) then begin
	jul = dblarr(abs(nday))
	sign = nday/abs(nday)
	jul = fjul + sign*dindgen(abs(nday))
	
	caldat,jul,mm,dd,yyyy,hh,mn,ss
	if (fflag eq 1) then $
		days = string(yyyy,format='(i4.4)')+string(mm,format='(i2.2)')+ $
		   	   string(dd,format='(i2.2)') $
	else $
		days = string(yyyy,format='(i4.4)')+'-'+string(mm,format='(i2.2)')+ $
		   	   '-'+string(dd,format='(i2.2)')
	
	last_day = days(abs(nday)-1)
	julian = jul
	return,days
endif else begin
	;Test of last_day as an input paramater
	if (~keyword_set(last_day)) then begin
		message,/INFO,'The parameter LAST_DAY must be specified!'
		return,''
	endif
	lday = strsplit(strtrim(last_day[0],2),'-',/EXTRACT)
	lflag = n_elements(lday)
	case lflag of
		1:begin
			yyyy = strmid(lday,0,4) & mm = strmid(lday,4,2) & dd = strmid(lday,6,2)
		end
		3:begin
			yyyy = lday(0) & mm = lday(1) & dd = lday(2)
		end
		else:begin
			message,/CONT,'Wrong input format for last_day!'
			return,''
		end
	endcase
	
	;Get the julian day for lday
    ljul = (julday(mm,dd,yyyy,hh,mn,ss))[0]
    
    ;Calulating the number of day(s)
    nday = round(ljul-fjul+1D)
    sign = nday/abs(nday)
	jul = fjul + sign*dindgen(abs(nday))
	
    caldat,jul,mm,dd,yyyy,hh,mn,ss
  
	if (fflag eq 1) then $
		days = string(yyyy,format='(i4.4)')+string(mm,format='(i2.2)')+ $
		   	   string(dd,format='(i2.2)') $
	else $
		days = string(yyyy,format='(i4.4)')+'-'+string(mm,format='(i2.2)')+ $
		   	   '-'+string(dd,format='(i2.2)')

    julian = jul
    return,days
endelse

END