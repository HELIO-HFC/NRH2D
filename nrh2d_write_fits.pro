;+
; NAME:
;		nrh2d_write_fits
;
; PURPOSE:
; 		Produces FITS format file containing the
;		NRH 2D image + results of source detections. 
;
; CATEGORY:
;		I/O
;
; GROUP:
;		NRHD2
;
; CALLING SEQUENCE:
;		IDL>nrh2d_write_fits, image, index, feat_str
;
; INPUTS:
;		image 	   - 2d array containing NRH image.
;		index      - a structure containing FITS image index.
;		feat_str   - a structure containing the feature parameters.		
;	
; OPTIONAL INPUTS:
;		write_fits - Define type of fits file to produce:
;						write_fits = 0 -> No fits file produced.
;						write_fits = 1 -> Only image is returned.
;						write_fits = 2 -> image with detection results only.
;										  (Default value).
;       outfname   - Scalar of string type containing the name of the output fit file.
;		output_dir - Scalar of string type specifying the output directory.
;
; KEYWORD PARAMETERS:
;		/SILENT - Quiet mode.
;
; OUTPUTS:
;		None.
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
;		cvtime
;		mwritefits
;
; EXAMPLE:
;		None.		
;
; MODIFICATION HISTORY:
;		Written by X.Bonnin, 28-SEP-2011.
;
;       06-JAN-2012, X.Bonnin:  Added outfname optional input.
;                               Renamed outputs to feat_str.									
;-

PRO nrh2d_write_fits,image,index,feat_str,$
                     outfname=outfname, $
					 output_dir=output_dir,$
					 write_fits=write_fits,$
					 error=error,$
					 SILENT=SILENT

;[1]:Initialization of the program
;[1]:=============================
;On_error,2

error = 1

;Check input parameters
if (n_params() lt 3) then begin
	message,/INFO,'Call is:'
	print,'nrh2d_write_fits,image,index,feat_str,$'
	print,'                 outfname=outfname, $'
	print,'                 output_dir=output_dir,$'
	print,'                 write_fits=write_fits,$'
	print,'                 error=error,/SILENT'
	return
endif

if (~keyword_set(write_fits)) then write_fits = 2
if (write_fits eq 0) then begin
	error = 0
	return
endif
SILENT = keyword_set(SILENT)

outflag = keyword_set(outfname)

if (~keyword_set(output_dir)) then cd,current=outdir else outdir = strtrim(output_dir[0],2)
;[1]:=============================


;[2]:Preparing image
;[2]:===============
if (not outflag) then begin
    date_obs = index.date_obs
    date = strsplit(date_obs,'T',/EXTRACT)
    time = strjoin((strsplit(date[1],':.',/EXTRACT))[0:2])
    date = strjoin(strsplit(date[0],'-',/EXTRACT))
endif

;flux -> byte
array = bytscl(image,min=min(image,/NAN),max=max(image,/NAN),top=255)

;mask = bytarr(nt,nf)
if (write_fits eq 2) then begin
	for i=0L,n_elements(feat_str.cc)-1L do begin
		if (strtrim(feat_str(i).cc,2) ne '') or (strtrim(feat_str(i).cc,2) ne 'NULL') then continue
		cc = feat_cc_extract(feat_str(i).cc,[feat_str(i).cc_x_pix,feat_str(i).cc_y_pix])
		for k=0L,n_elements(cc[0,*])-1L do array[cc[0,k],cc[1,k]] = 255 ;array[cc[0,k],cc[1,k]]
	endfor
	;array = mask
	if (not outflag) then outfname = 'nrh2d_'+date+'T'+time+'_products.fits'
endif else begin
    if (not outflag) then outfname = 'nrh2d_'+date+'T'+time+'_image.fits'
endelse
;[2]:===============

;[3]:Preparing fits header
;[3]:=====================
outpath = outdir + path_sep() + outfname
mwritefits,index,array,outfile=outpath
;[3]:=====================

error = 0
END
