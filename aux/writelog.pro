;+
; NAME:
;		writelog
;
; PURPOSE:
; 		Writes a message in a LOG file.
;
; CATEGORY:
;		I/O
;
; GROUP:
;		None.
;
; CALLING SEQUENCE:
;		IDL> writelog,lun,msg
;
; INPUTS:
;		lun - Logical unit number of the LOG file.
;		msg	- Message to write in the LOG file (scalar of string type.)	
;	
; OPTIONAL INPUTS:
;		filename - Scalar of string type containing the path name of the LOG file to create.
;				   (filename must be set to create the LOG file.)
;
; KEYWORD PARAMETERS:
;		/OPEN	- open log file.
;		/CLOSE	- close log file.
;		/SILENT	- quiet mode : the message is only written in the LOG file 
;				  (but not printed on the terminal screen.)
;		/NOLOG	- Do not write log file (if is set, the routine will just print information on prompt).
;		/HELP 	- Display help.
;
; OUTPUTS:
;		None.
; OPTIONAL OUTPUTS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; CALL:
;		None.
; EXAMPLE:
;		None.	
;
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin, 28-OCT-2010. 
;
;		X.Bonnin, 03-JUN-2011:	Added /OPEN and /CLOSE keywords.
; 
;-

PRO writelog,lun,msg, $
		     filename=filename, $
			 OPEN=OPEN,CLOSE=CLOSE,$
			 SILENT=SILENT,HELP=HELP


;[1]: Initializing the program
;==============================
if (keyword_set(HELP)) then begin
	message,/INFO,'Call is:'
	print,'writelog,lun,msg,filename=filename,/OPEN,/CLOSE,/SILENT,/HELP'
	return
endif

file_set = keyword_set(filename)
lun_set = keyword_set(lun)
msg_set = n_elements(msg)
SILENT = keyword_set(SILENT)
OPEN = keyword_set(OPEN)
CLOSE = keyword_set(CLOSE)
;==============================

;[2]: Editing LOG file
;=====================

;Openning LOG file----
if (file_set) and (OPEN) then begin

	;Defining name of the log file (TBD)
	logfilename = strtrim(filename[0],2)
	
	;Openning log file
	openw,lun,logfilename,error=err,/GET_LUN
	if (err ne 0) then message,!error_state.msg
	return
endif else lun = 0

;Closing LOG file-----
if (lun_set) and (CLOSE) then begin
	close,lun
	free_lun,lun
	return
endif

;Writting in LOG file-
if (msg_set) then begin
	mess = strtrim(msg[0],2)
	
	if (lun_set) then printf,lun,mess
	if (~SILENT) then print,mess
	
	return
endif
;=====================


END
