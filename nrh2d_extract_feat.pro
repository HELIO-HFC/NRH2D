;+
; NAME:
;		nrh2d_extract_feat
;
; PURPOSE:
; 		Extract Gaussian radio sources from 2D NRH image.
;
; CATEGORY:
;		Image processing
;
; GROUP:
;		NRH2D
;
; CALLING SEQUENCE:
;		feat_str_out = nrh2d_extract_feat(image,index,feat_str_in)
;
; INPUTS:
;		image	    - 2D array containing the NRH image to process.
;		index       - Structure containing the corresponding fits header.
;		feat_str_in - Structure to fill with feature parameters extracted.
;                     (See hfc_radiosources__define.pro routine for
;                     more details about fields in structure). 
;	
; OPTIONAL INPUTS:
;		fact  - See extract_sources.pro header. Default is 1.
;       nbmax - See extract_sources.pro header. Default is 3.
;       seuil - See extract_sources.pro header. Default is 2.
;       fmax  - See extract_sources.pro header. Default is 3.
;       amax  - Maximum area (in percentage of the full solar disk area) authorized
;               for the detected features. 
;               (If a feature's area is greater than amax, it will be ignored.)
;               Default is 0.5
;
; KEYWORD PARAMETERS:
;		/SILENT      - Quiet mode.
;
; OUTPUTS:
;		feat_str_out - feat_str_in updated with feature parameters.		
;
; OPTIONAL OUTPUTS:
;		features  - Feature parameters extracted by the algorithm.
;		error     - Return 1 if an error occurs, 0 else.
;		
; COMMON BLOCKS:		
;		None.
;	
; SjE EFFECTS:
;		None.
;		
; RESTRICTIONS/COMMENTS:
;		The SSW library including the NRH package must be
;		loaded.
;			
; CALL:
;		extract_sources
;       feat_cc_make
;       feat_cc_extract
;		hfc_pix2car
;
; EXAMPLE:
;		None.		
;
; MODIFICATION HISTORY:
;		Written by X.Bonnin.			
;				
;-

FUNCTION nrh2d_extract_feat,image,index,feat_str_in,$
                            nbmax=nbmax, fact=fact, seuil=seuil, $
                            fmax=fmax, amax=amax, $
			                features=features,error=error,$
		                    SILENT=SILENT

error = 1
;[1]:Initialize the program
;[1]:======================
if (n_params() lt 2) then begin
	message,/INFO,'Call is:'
	print,'feat_str_out = nrh2d_extract_feat(image,index,feat_str_in,$'
    print,'                                  nbmax=nbmax, fact=fact, seuil=seuil, $'
    print,'                                  fmax=fmax, amax=amax, $ '
	print,'                                  features=features,$'
	print,'                                  error=error,/SILENT)'
	return,0
endif

if (not keyword_set(nbmax)) then nbmax = 3
if (not keyword_set(fact)) then fact = 1
if (not keyword_set(seuil)) then seuil = 2
if (not keyword_set(fmax)) then fmax = 0.8
if (not keyword_set(amax)) then amax = 0.5

;Convert degrees to megameter on the Sun surface
deg2mm = 6.96e8/!radeg 

features = {gra:0.,grb:0.,tet:0.,max:0.,xmax:0.,ymax:0.}
features = replicate(features,nbmax)

feat_str_out = feat_str_in

SILENT = keyword_set(SILENT)
;[1]:======================

;[2]:Extract radio sources
;[2]:=====================

extract_sources, image, features, FACT=fact, NBMAX=nbmax, SEUIL=seuil, FMAX=fmax
test = (features(*).max NE -1.0) AND (features(*).max NE 0.0)
goodvalues = where(test EQ 1, count)

if (count NE 0) then begin
	feat_str_out = replicate(feat_str_out,count)

	;display2d,image,color_table=3
    ;loadct,39,/SILENT
 	;col = bytscl(lindgen(count),top=150)+20
    
    for j=0L, count-1L do begin
      	;oplot,features(j).xmax+fltarr(2),features(j).ymax + fltarr(2),psym=1,color=col(j),symsize=2.,thick=2.
      	
      	;Store ellipse parameters
      	feat_str_out(j).el_axis1 = sqrt(1./features(j).gra) ;long axis (in pix)
      	feat_str_out(j).el_axis2 = sqrt(1./features(j).grb) ;short axis (in pix))
      	feat_str_out(j).el_angle = features(j).tet*!radeg   ;direction angle (in deg)
    
      	;Store gravity center coordinates
      	feat_str_out(j).feat_x_pix = features(j).xmax
      	feat_str_out(j).feat_y_pix = features(j).ymax
      	feat_str_out(j).feat_x_arcsec = index.cdelt1*(features(j).xmax - index.crpix1)
      	feat_str_out(j).feat_y_arcsec = index.cdelt2*(features(j).ymax - index.crpix2)
      	Xcar = hfc_pix2car([feat_str_out(j).feat_x_pix,feat_str_out(j).feat_y_pix],$
      						index.cdelt1, index.cdelt2,$
      						index.crpix1,index.crpix2,$
      						index.solar_r,index.date_obs,$
      						hel_coord=Xhel)
      	feat_str_out(j).feat_hg_long_deg = Xhel[0]
      	feat_str_out(j).feat_hg_lat_deg = Xhel[1]
      	feat_str_out(j).feat_carr_long_deg = Xcar[0]
      	feat_str_out(j).feat_carr_lat_deg = Xcar[1]
      	
      	
      	;Compute radio source contours
      	theta = 2.*!pi*findgen(361)/360.
		xell = sqrt(1./features(j).gra)*cos(theta)
		yell = sqrt(1./features(j).grb)*sin(theta)
		x = features(j).xmax + xell*cos(features(j).tet) - yell*sin(features(j).tet) ;pix
		y = features(j).ymax + xell*sin(features(j).tet) + yell*cos(features(j).tet) ;pix
      	x = round(x)
      	y = round(y)
      	
      	;oplot,x,y,col=col(j),line=2
 		
 		;Get feature contour
 		pix_inside = polyfillv(x,y,index.naxis1, index.naxis2)
		mask = intarr(index.naxis1, index.naxis2)
		mask[pix_inside] = 1
		kernel = intarr(3,3) + 1
		mask = mask - erode(mask,kernel)
	
		;Generate chain code of contour
		cc = feat_cc_make(mask,start_pix=start_pix,error=error)

		feat_str_out(j).cc = cc
		feat_str_out(j).cc_length = strlen(cc)
		feat_str_out(j).cc_x_pix = start_pix[0]
		feat_str_out(j).cc_y_pix = start_pix[1]
		feat_str_out(j).cc_x_arcsec = index.cdelt1*(start_pix[0] - index.crpix1)
		feat_str_out(j).cc_y_arcsec = index.cdelt2*(start_pix[1] - index.crpix2)
		
		;Generate chain code boundaries
		bd = feat_cc_extract(cc,start_pix)
		nbd = n_elements(bd[0,*])
		
		;oplot,bd[0,*],bd[1,*],line=1,col=50,thick=2.5
	
		;Get bounding rectangle corners coordinates
		Xmin = min(bd[0,*],max=Xmax)
		Ymin = min(bd[1,*],max=Ymax)
		
		;in pixel and arcsec
		feat_str_out(j).br_x0_pix=Xmin
		feat_str_out(j).br_y0_pix=Ymin
	 	feat_str_out(j).br_x1_pix=Xmin
	 	feat_str_out(j).br_y1_pix=Ymax
	 	feat_str_out(j).br_x2_pix=Xmax
	 	feat_str_out(j).br_y2_pix=Ymin
	 	feat_str_out(j).br_x3_pix=Xmax
	 	feat_str_out(j).br_y3_pix=Ymax
	 	feat_str_out(j).br_x0_arcsec=index.cdelt1*(Xmin - index.crpix1)
	 	feat_str_out(j).br_y0_arcsec=index.cdelt2*(Ymin - index.crpix2)
	 	feat_str_out(j).br_x1_arcsec=index.cdelt1*(Xmin - index.crpix1)
	 	feat_str_out(j).br_y1_arcsec=index.cdelt2*(Ymax - index.crpix2)
	 	feat_str_out(j).br_x2_arcsec=index.cdelt1*(Xmax - index.crpix1)
	 	feat_str_out(j).br_y2_arcsec=index.cdelt2*(Ymin - index.crpix2)
	 	feat_str_out(j).br_x3_arcsec=index.cdelt1*(Xmax - index.crpix1)
	 	feat_str_out(j).br_y3_arcsec=index.cdelt2*(Ymax - index.crpix2)
	 	
	 	;oplot,x,y,color=col(j)
		;oplot,[Xmin,Xmin],[Ymin,Ymax],color=col(j)
		;oplot,[Xmax,Xmax],[Ymin,Ymax],color=col(j)
		;oplot,[Xmin,Xmax],[Ymin,Ymin],color=col(j)
		;oplot,[Xmin,Xmax],[Ymax,Ymax],color=col(j)

	 	;In heliographic and carrington (degrees)
	 	Xc = fltarr(2,nbd) & Xh = fltarr(2,nbd)
	 	for k=0L,nbd-1L do begin
	 		Xc_k = hfc_pix2car(reform(bd[*,k]),$
	      						index.cdelt1, index.cdelt2,$
	      						index.crpix1,index.crpix2,$
	      					  	index.solar_r,index.date_obs,$
	      					  	hel_coord=Xh_k)
	      					  
	      	Xc[*,k] = Xc_k
	      	Xh[*,k] = Xh_k
		endfor
		
		Xhmin = min(Xh[0,*],max=Xhmax)
		Xcmin = min(Xc[0,*],max=Xcmax)
		Ymin = min(Xc[1,*],max=Ymax)
		
	 	feat_str_out(j).br_hg_long0_deg=Xhmin
	 	feat_str_out(j).br_hg_lat0_deg=Ymin
	 	feat_str_out(j).br_hg_long1_deg=Xhmin
	 	feat_str_out(j).br_hg_lat1_deg=Ymax
	 	feat_str_out(j).br_hg_long2_deg=Xhmax
	 	feat_str_out(j).br_hg_lat2_deg=Ymin
	 	feat_str_out(j).br_hg_long3_deg=Xhmax
	 	feat_str_out(j).br_hg_lat3_deg=Ymax
	 	feat_str_out(j).br_carr_long0_deg=Xcmin
	 	feat_str_out(j).br_carr_lat0_deg=Ymin
	 	feat_str_out(j).br_carr_long1_deg=Xcmin
	 	feat_str_out(j).br_carr_lat1_deg=Ymax
	 	feat_str_out(j).br_carr_long2_deg=Xcmax
	 	feat_str_out(j).br_carr_lat2_deg=Ymin
	 	feat_str_out(j).br_carr_long3_deg=Xcmax
	 	feat_str_out(j).br_carr_lat3_deg=Ymax


      	;Compute max flux (sfu)
      	feat_str_out(j).feat_max_int = features(j).max
      	;Compute mean flux (sfu)
      	feat_str_out(j).feat_mean_int = mean(image[pix_inside])
      	
      	;Compute area
      	feat_str_out(j).feat_area_pix = n_elements(pix_inside)
      	feat_str_out(j).feat_area_deg2 = poly_area(bd[0,*],bd[1,*])
      	feat_str_out(j).feat_area_mm2  = feat_str_out(j).feat_area_deg2*(deg2mm)^2
      	
      	A_sun = long(!pi*index.solar_r^2)
 ;     	print,inputs.amax*A_sun,feat_str_out(j).feat_area_pix
 
      	if (feat_str_out(j).feat_area_pix gt amax*A_sun) then continue
    	feat_str_out(j).id_rs = j+1L

  	endfor
endif
;[2]:=====================

error = 0
return,feat_str_out
END