;+
; NAME:
;		nrh2d_hfc_struct
;
; PURPOSE:
; 		Load fields values of HFC tables into structures.
;
; CATEGORY:
;		I/O
;
; GROUP:
;		RABAT3
;
; CALLING SEQUENCE:
;		IDL> nrh2d_hfc_struct,str_obs,str_oby,str_frc,str_feat
;
; INPUTS:
;		None.
;	
; OPTIONAL INPUTS:
;		None.
;
; KEYWORD PARAMETERS:
;		/SILENT - Quiet mode.
;
; OUTPUTS:
;		str_obs  - Empty structure containing the fields of the OBSERVATIONS HFC table. 
;		str_oby  - Empty structure containing the fields of the OBSERVATORY HFC table.
;		str_frc  - Empty structure containing the fields of the FRC_INFO HFC table.
;		str_feat - Empty structure containing the fields of the RADIOSOURCES HFC table.
;
; OPTIONAL OUTPUTS:
;		error - Equal to 1 if an error occurs, 0 else.
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
;		hfc_observations__define
;		hfc_observatory__define
;		hfc_frc_info__define
;		hfc_radiosources__define		
;
; EXAMPLE:
;		None.
;
; MODIFICATION HISTORY:
;		Written by:		X.Bonnin, 06-JAN-2012.
;
;-

PRO nrh2d_hfc_struct,str_obs,str_oby,str_frc,str_feat, $
					  error=error
						   
	
;Intializing structures
str_obs = {hfc_observations}
str_oby = {hfc_observatory}
str_frc = {hfc_frc_info}
str_feat = {hfc_radiosources}

error = 0
return
END
