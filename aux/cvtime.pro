;+
; NAME:
;		cvtime
;
; PURPOSE:
; 		Converts time in float format to string 
;		(or string to float).
;
; CATEGORY:
;		I/O
;
; GROUP:
;		RABAT
;
; CALLING SEQUENCE:
;		IDL>Results = cvtime(time)
;
; INPUTS:
;		time - Scalar or vector containing the time(s) to convert.
;			   time can be of float, double or string type.
;			   if time is of string type, input format must be 'hh:mm:ss'.	
;	
; OPTIONAL INPUTS:
;		None.
;
; KEYWORD PARAMETERS:
;		/DOUBLE	- Use double type instead of float.
;
; OUTPUTS:
;		new_time - Input time converted in float type (double if keyword /DOUBLE is set), 
;				   if the input time is of string type.
;                  Input time converted in string type (format 'hh:mm:ss'),
;				   if the input time is of float/double type.
;		
; OPTIONAL OUTPUTS:
;		None.
;		
; COMMON BLOCKS:		
;		None.	
;	
; SIDE EFFECTS:
;		None.
;		
; RESTRICTIONS/COMMENTS:
;		None.
;		
; CALL:
;		None.
;
; EXAMPLE:
;		None.		
;
; MODIFICATION HISTORY:
;		Written by X.Bonnin, 05-APR-2011.
;									
;-

FUNCTION cvtime,time,DOUBLE=DOUBLE

if (n_params() lt 1) then begin
	message,/INFO,'Call is:'
	print,'Results = cvtime(time,/DOUBLE)'
	return,0
endif

DOUBLE = keyword_set(DOUBLE)
n = n_elements(time)

type = size(time,/TYPE)
if (type ne 7) then begin
	hh = fix(time)
	mm = 60.0d*(time - hh)
	ss = 60.0d*(mm - fix(mm))

	new_time = string(hh,format='(i2.2)')+':'+$
		   		string(fix(mm),format='(i2.2)')+':'+$
		   		string(fix(ss),format='(i2.2)')
endif else begin
	new_time = dblarr(n)
	for i=0L,n-1L do begin
		new_time_i = double(strsplit(time[i],':',/EXTRACT))
		new_time[i] = double(new_time_i[0]) + $
					  double(new_time_i[1])/60.d0 + $
					  double(new_time_i[2])/3600.d0
	endfor
	if (~DOUBLE) then new_time = float(new_time)
endelse

return,new_time
END