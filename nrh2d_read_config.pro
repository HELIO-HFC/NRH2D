FUNCTION nrh2d_read_config,config_file, $
                           error=error, $
					       SILENT=SILENT

;+
; NAME:
;		nrh2d_read_config
;
; PURPOSE:
; 		Read input parameters written in the
;		configuration file, and update info
;		structures.
;
; CATEGORY:
;		I/O
;
; GROUP:
;		NRH2D
;
; CALLING SEQUENCE:
;		IDL> inputs_str = nrh2d_read_config(config_file)
;
; INPUTS:
;		config_file - Full path name to the configuration file 
;                     containing the input parameters (scalar of string type).
;	
; OPTIONAL INPUTS:
;		None.
;
; KEYWORD PARAMETERS:
;		/SILENT	- Quiet mode.
;
; OUTPUTS:
;		inputs_str - Structure containing the input parameters for NRH2D.
;                    See nrh2d_inputs__define.pro header for the list of
;                    fields returned. 
;
; OPTIONAL OUTPUTS:
;		error - Equal to 1 if an error occurs, 0 otherwise.
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
;		nrh2d_config__define		
;
; EXAMPLE:
;		None.
;
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin, 26-SEP-2011.
;
;       15-DEC-2011, X.Bonnin:  Added nrh2_config__define call.
;-

error = 1
if (n_params() lt 1) then begin
	message,/INFO,'Call is:'
	print,'inputs_str = nrh2d_read_config(config_file ,$'
	print,'                               ,error=error,/SILENT)'
	return,''
endif

SILENT = keyword_set(SILENT)

in_file = strtrim(config_file[0],2)
if (~file_test(in_file)) then message,'No input file found!'

inputs_str = {nrh2d_config}
ntags = n_tags(inputs_str)
nlines  = file_lines(in_file)

if (ntags ne nlines) then message,'Incompatible numbers of fields!'

tags = strarr(nlines) 
data = strarr(nlines)
openr,lun,in_file,/GET_LUN
for i=0,nlines-1 do begin
	data_i = ""
	readf,lun,data_i
	data_i = strtrim(strsplit(data_i,'=',/EXTRACT),2)
	tags[i] = strlowcase(data_i[0])
	data[i] = data_i[1]
	jflag = execute('inputs_str.'+tags[i]+'=+data[i]')
endfor
close,lun
free_lun,lun

error = 0
return,inputs_str
END
