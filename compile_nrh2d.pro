;+
; NAME:
;		compile_nrh2d
;
; PURPOSE:
; 		Compiles all the nrh2d IDL routines.
;
; CATEGORY:
;		IDL batch file 
;		
; GROUP:
;		NRH2D
;
; CALLING SEQUENCE:
;		IDL> @compile_nrh2d
;
; INPUTS:
;		None. 
;	
; OPTIONAL INPUTS:
;		None.	
;
; KEYWORD PARAMETERS:
;		None.
;
; OUTPUTS:
;		None.		
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
; RESTRICTIONS:
;		None.
;
; CALL:
;		None.		
;
; EXAMPLE:
;		IDL>@compile_nrh2d
;
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin, 26-SEP-2011.
;
;-

.compile nrh2d_hfc_struct
.compile nrh2d_config__define
.compile nrh2d_read_config
.compile nrh2d_extract_feat
.compile nrh2d_write_fits
.compile nrh2d
