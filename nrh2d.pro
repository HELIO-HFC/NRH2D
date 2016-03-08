PRO nrh2d,data_file,config_file,$
		  freq=freq,cadence=cadence,$
		  hbeg=hbeg,hend=hend,$
		  data_dir=data_dir,$
		  output_dir=output_dir,$
		  input_dir=input_dir,$
          frc_struct=frc_stc, oby_struct = oby_stc, $
          obs_struct=obs_struct, feat_struct=feat_struct, $
          GET_CLOSEST_FREQ=GET_CLOSEST_FREQ, $ 
		  write_png=write_png, WRITE_FITS=WRITE_FITS,$
		  WRITE_CSV=WRITE_CSV,$
		  WRITE_LOG=WRITE_LOG,$
		  SILENT=SILENT


;+
; NAME:
;		nrh2d
;
; PURPOSE:
; 		This procedure detects gaussian 
;		radio sources on a 
;       Nancay RadioHeliograph (NRH) 2D image.
;       To perform extraction, it uses dedicated routines of the 
;		Solar SoftWare (SSW) NRH package.
;
; CATEGORY:
;		Image processing
;
; GROUP:
;		NRH2D
;
; CALLING SEQUENCE:
;		IDL>nrh2d,data_file,config_file
;
; INPUTS:
;		data_file   - Scalar of string type specifying the name of
;				      input NRH data file to process.
;				      (The procedure reads fits format file only.)
;		config_file - Scalar of string type 
;					  providing the name of the configuration file
;                     where input parameters are written.  
;	
; OPTIONAL INPUTS:
;		freq		- Frequency of observations:
; 							freq = 228 MHz
; 							freq = 298 MHz
;							freq = 360 MHz
; 							freq = 327 MHz
; 							freq = 270 MHz
; 							freq = 445 MHz
; 							freq = 408 MHz
; 							freq = 150 MHz
; 							freq = 173 MHz
; 							freq = 432 MHz
;                           freq = 164 MHz
;					  Defaults is 150 MHz.
;		hbeg        - Scalar of string type containg the time of the first image to read (in UTC).
;					  Default is the first image available.
;		hend    	- Scalar of string type containg the time of the last image to read. 
;					  Default is '15:00:00'.
;		data_dir    - Scalar of string type providing 
;				      the path to the data file directory.
;		output_dir  - Scalar of string type containing the path to the 
;					  directory where output data files will be saved.
;					  (The default is current directory.)		
;		input_dir 	- Scalar of string type providing 
;					  the path to the inputs file directory.
;		cadence	    - Cadence (in seconds) of the images processed.
;					  Default is 10 seconds (highest resolution).
;		write_png   - Write output png image files:
;							write_png = 0 --> No output
;							write_png = 1 --> Write png 2d image of the radio observation
;											  without detection results.
;							write_png = 2 --> Write png 2d imageof the radio obervation
;											  including detection results.
;
; KEYWORD PARAMETERS:
;		/WRITE_CVS        - Write cvs files containing detection products.
;		/WRITE_FITS	      - Write fits file of the 2d image with detection
;					        results overplotted.
;       /GET_CLOSEST_FREQ - If the input frequency is not available, get the closet one available. 
;		/SILENT           - Quiet mode.
;
; OUTPUTS:
;		None.		
;
; OPTIONAL OUTPUTS:
;		frc_struct  - Structure containing frc info.
;       oby_struct  - Structure containing observatory info.
;       obs_struct  - Structure containing observations info.
;       feat_struct - Structure containing feature parameters extracted.
;		
; COMMON BLOCKS:		
;		CLOG	
;	
; SIDE EFFECTS:
;		None.
;		
; RESTRICTIONS/COMMENTS:
;		The SSW library with the NRH package must be
;		loaded.
;			
; CALL:
;		anytim
;		read_nrh
;		cvtime
;		tim2carr
;		jd2str
;		anytim2jd
;       hfc_write_csv
;       nrh2d_hfc_struct
;       nrh2d_read_config
;		nrh2d_extract_feat
;		nrh2d_write_csv
;		nrh2d_write_fits
;		write_posinrh
;       rh_open
;       rh_close
;
; EXAMPLE:
;		None.		
;
; MODIFICATION HISTORY:
;		Written by X.Bonnin.			
;				
; Version 1.00
;		29-SEP-2011, X.Bonnin:	First release.
;
; Version 1.01
;		28-NOV-2011, X.Bonnin:  Added cadence, hbeg, and hend optional inputs.
;
; Version 1.02
;		20-DEC-2011, X.Bonnin:	Fixed a minor bug concerning the output format.
;								Changed the nomenclature of feature output files (uses feat instead of par).
; Version 1.03
;       06-JAN-2012, X.Bonnin:  Introduced the use of hfc structures definition.
;                               Renamed config_file to inputs_file.
;                               Renamed nrh2d_readconf routine to nrh2d_read_config.
;                               Renamed nrh2d_extract routine to nrh2d_extract_feat.
;                               Renamed nrh2d_writecsv/fits routine to nrh2d_write_csv/fits.
; Version 1.04
;		14-JAN-2012, X.Bonnin:	Added write_png optional input.
; Version 1.05
;       26-MAR-2012, X.Bonnin:  Load frequency channel from the input file using rh_open.pro 
;                               to check input frequency avaibility.
;                               Added /GET_CLOSEST_FREQ optional keyword.
;                               Radio flux is provided in SFU (10^-22 W/M^2/Hz).
;
;-
version = '1.05'

;[1]:Initialize the program
;[1]:======================
if (n_params() lt 2) then begin
	message,/INFO,'Call is:'
	print,'nrh2d,data_file,config_file,$'
	print,'      freq=freq,cadence=cadence,$'
	print,'      hbeg=hbeg,hend=hend,$'
	print,'      data_dir=data_dir,$'
	print,'      output_dir=output_dir,$'
	print,'      input_dir=input_dir,$'
    print,'      frc_struct=frc_struct, oby_struct=oby_struct, $'
    print,'      obs_struct=obs_struct, feat_struct=feat_struct, $'
	print,'      write_png=write_png, /WRITE_FITS,$'
	print,'      /WRITE_CSV,/WRITE_LOG,$'
	print,'      /GET_CLOSEST_FREQ,/SILENT'
	return
endif
syst0 = systime(1)

WRITE_FITS = keyword_set(WRITE_FITS)
WRITE_CSV= keyword_set(WRITE_CSV)
WRITE_LOG = keyword_set(WRITE_LOG)
GET_CLOSEST_FREQ = keyword_set(GET_CLOSEST_FREQ)
SILENT = keyword_set(SILENT)

if (not keyword_set(write_png)) then wpng = 0 else begin
	wpng = fix(write_png[0])
	loadct,3,/SILENT ;red temp color table
	tvlct,r,g,b,/GET
endelse

;constant to convert julian days to seconds
jd2sec = 24.0d*3600.0d

cd,current=current_dir

dfile = strtrim(data_file[0],2)
if (keyword_set(data_dir)) then dfile = strtrim(data_dir[0],2)+path_sep()+file_basename(dfile)
dfile = (file_search(dfile))[0]
if (dfile eq '') then message,'Data file '+dfile+' not found!'
data_dir = file_dirname(dfile)

ifile = strtrim(config_file[0],2)
if (keyword_set(input_dir)) then ifile = strtrim(input_dir[0],2)+path_sep()+file_basename(ifile)
ifile = (file_search(ifile))[0]
if (ifile eq '') then message,'Configuration file '+ifile+' not found!'
input_dir = file_dirname(ifile)

if (~keyword_set(output_dir)) then begin
	output_dir = current_dir + path_sep() + 'Products'
	if (~file_test(output_dir,/DIR)) then spawn,'mkdir '+output_dir	
endif

outfnroot = 'nrh2d_'+strjoin(strsplit(version,'.',/EXTRACT))

;Hostname
spawn,'hostname',host

;Current date
run_date = (strsplit(anytim(!stime, /ccsds),'.',/EXTRACT))[0]

;Create LOG file (if WRITE_LOG is set)
logname = output_dir+path_sep()+outfnroot+strjoin(strsplit(run_date,'-T:.',/EXTRACT))+'.log'
writelog,lunlog,filename=logname,OPEN=WRITE_LOG

;Writting LOG file header
writelog,lunlog,'Starting NRH2D software '+version
writelog,lunlog,'run_date :'+run_date
writelog,lunlog,'host: '+host

;Cadence in sec
max_cad = 10.0d
if (~keyword_set(cadence)) then cad = max_cad else cad = double(cadence[0])

if (~keyword_set(freq)) then f = '150' else f = strtrim(freq[0],2)

;Get frequency list
@rh_common.inc
if not (rh_open(dfile)) then begin
    message,/CONT,'Cannot read input data file!'
    return
endif
rh_close

wave = f
freqList = entfi.frq/10
channel = (where(fix(f) eq freqList))[0]
if (channel eq -1) then begin
    message,/CONT,'Input Frequency ('+f+') not available!'
    if (GET_CLOSEST_FREQ) then begin
        ch_min = min(abs(fix(f) - freqList),channel)
        f = strtrim(freqList[channel],2)
        wave = [wave,f]
        print,'Get the closest frequency channel available: '+f
    endif else begin    
        print,entfi.frq/10
        return
    endelse
endif

if (~keyword_set(hbeg)) then hbeg = strjoin(string(entfi.hdeb,format='(i2.2)'),':')
if (~keyword_set(hend)) then hend = strjoin(string(entfi.hfin,format='(i2.2)'),':')

writelog,lunlog,'date to process: '+strjoin([string(entfi.dat[0:1],format='(i2.2)'),strtrim(entfi.dat[2],2)],'-')
writelog,lunlog,'Start time: '+hbeg
writelog,lunlog,'End time: '+hend
writelog,lunlog,'Channel: '+f+' ('+strtrim(channel,2)+')'

;Initializing structure 
writelog,lunlog,'--> Loading input parameters...'
in_stc = nrh2d_read_config(ifile)
if (size(in_stc,/TNAME) ne 'STRUCT') then begin
    writelog,lunlog,'ERROR: cannot read input parameters correctly!'
    return
endif
if (strtrim(version,2) ne strtrim(in_stc.version,2)) then writelog,lunlog,'WARNING:Incompatible software version!' 

nrh2d_hfc_struct,obs_stc,oby_stc,frc_stc,feat_stc

;Update frc info
frc_stc.id_frc_info = 1
frc_stc.institut = in_stc.institut
frc_stc.code = in_stc.code
frc_stc.version = version
frc_stc.feature_name = in_stc.feature_name
frc_stc.enc_met = in_stc.enc_met
frc_stc.person = in_stc.person
frc_stc.contact = in_stc.contact
frc_stc.reference = in_stc.reference

;Update observatory info
nwave = n_elements(wave)
oby_stc = replicate(oby_stc,nwave)
oby_stc.id_observatory = 1 + indgen(nwave)
oby_stc.observat = in_stc.observat + strarr(nwave)
oby_stc.instrume = in_stc.instrume + strarr(nwave)
oby_stc.telescop = in_stc.telescop + strarr(nwave)
oby_stc.units = 'SFU' + strarr(nwave)
oby_stc.wavemin = wave
oby_stc.wavemax = wave
oby_stc.wavename = in_stc.wavename  + strarr(nwave)
oby_stc.waveunit = 'MHz'  + strarr(nwave)
oby_stc.spectral_name = in_stc.spectral_name  + strarr(nwave)
oby_stc.obs_type = in_stc.obs_type  + strarr(nwave)

;Update observations info
obs_stc.id_observations = 1
obs_stc.observatory_id = max(oby_stc.id_observatory)
obs_stc.filename = file_basename(dfile)
obs_stc.file_format = "Binary"
obs_stc.loc_filename = data_dir + path_sep() + file_basename(dfile)
obs_stc.comment = 'NRH 1D to 2D interferometer deconvolution was performed using the SSW/NRH package.'

;Update feature struct
feat_stc.frc_info_id = 1
feat_stc.observations_id = 1
feat_stc.run_date = run_date

writelog,lunlog,'--> Loading input parameters:OK'
;[1]:======================

;[2]:Writing frc_info and observatory output files
;[2]:=============================================
if (WRITE_CSV) then begin
    writelog,lunlog,'--> Writing frc info output file...'
    frc_file = outfnroot+'_nan_frc_info.csv'
    frc_path = output_dir + path_sep() + frc_file
    hfc_write_csv,frc_stc,frc_path
    writelog,lunlog,'--> Writing frc info output file:OK'
    
    writelog,lunlog,'--> Writing observatory output file...'
    oby_file = outfnroot+'_nan_observatory.csv'
    oby_path = output_dir + path_sep() + oby_file
    nlines = 0
    if (file_test(oby_path,/REG)) then nlines = file_lines(oby_path)
    if (nwave gt nlines-1) then hfc_write_csv,oby_stc,oby_path else $
        print,oby_path+' already exists.'
    writelog,lunlog,'--> Writing observatory output file:OK'
endif
;[2]:=============================================

;[3]:Load nrh fits file
;[3]:==================
writelog,lunlog,'--> Reading '+dfile+'...'
read_nrh,file_basename(dfile),index,data,freq=channel,$
		 hbeg=hbeg,hend=hend,dir=data_dir,/SFU
nimg = n_elements(index)
writelog,lunlog,'Number of image(s) found: '+strtrim(nimg,2)
if (nimg eq 0) then begin
	writelog,lunlog,'ERROR:Empty nrh data file!'
	return
endif
date_obs = index.date_obs
jd_obs = anytim2jd(date_obs) & jd_obs = jd_obs.int + jd_obs.frac
jd_obs0 = min(jd_obs,imin,max=jd_obs1,subscript_max=imax,/NAN)
starttime =  date_obs[imin]
endtime = date_obs[imax]
jd_cad = cad/jd2sec
jd_cad_max = max_cad/jd2sec
nimg = long((jd_obs1 - jd_obs0)/jd_cad) + 1L
writelog,lunlog,strtrim(nimg,2)+' images to process between '+$
		 starttime+' and '+endtime
writelog,lunlog,'--> Reading '+dfile+':OK'
;[3]:==================

;[4]:Extract and save radio sources
;[4]:==============================
writelog,lunlog,'--> Running detection...'
iter = -1l
jd_obs_i = jd_obs0 & ni = 0L
obs_struct = obs_stc & feat_struct = feat_stc
while (jd_obs_i lt jd_obs1 + 0.5d*jd_cad_max) do begin
    iter++

	dt_i = min(abs(jd_obs_i - jd_obs),i)
	date_obs_i = jd2str(jd_obs_i,format=1)
    cdate_i = strjoin(strsplit(date_obs_i,':-',/EXTRACT))
	writelog,lunlog,strtrim(nimg-iter,2)+': Searching an image around '+date_obs_i
	if (dt_i gt jd_cad_max) then begin
		writelog,lunlog,'No image found!'
		jd_obs_i = jd_obs_i + jd_cad
		continue
	endif 
	ni++
	im_i = reform(data[*,*,i])
	ind_i = index(i)

	writelog,lunlog,'Image '+strtrim(i,2)+' found at '+date_obs[i]+' ['+strtrim(ni,2)+']'
	
    writelog,lunlog,'--> Updating observations meta-data...'
	obs_stc_i = obs_stc
    obs_stc_i.date_obs = date_obs_i
    obs_stc_i.date_end = jd2str(jd_obs_i + ind_i.exptime/jd2sec)
    jd = anytim2jd(date_obs_i) 
    obs_stc_i.jdint = jd.int
    obs_stc_i.jdfrac = jd.frac
    obs_stc_i.exp_time = index(i).exptime
    obs_stc_i.c_rotation = fix(tim2carr(date_obs_i,/DC))
    obs_stc_i.bitpix = ind_i.bitpix
    obs_stc_i.naxis1 = ind_i.naxis1
    obs_stc_i.naxis2 = ind_i.naxis2
    obs_stc_i.r_sun = ind_i.solar_r
    obs_stc_i.center_x = ind_i.crpix1
    obs_stc_i.center_y = ind_i.crpix2
    obs_stc_i.cdelt1 = ind_i.cdelt1
    obs_stc_i.cdelt2 = ind_i.cdelt2
    writelog,lunlog,'--> Updating observations meta-data:OK'    

	if (wpng gt 0) then begin
		print,'--> Writing 2d observation into a png image file...'
		imb_i = bytscl(im_i,top=255,/NAN)
		png_file = file_basename(obs_stc_i.filename)
		png_time = strtrim(strjoin((strsplit(obs_stc_i.date_obs,'T:',/EXTRACT))[1:*]),2)
		pos = strpos(png_file,'.',/REVERSE_S)
		png_file = strmid(png_file,0,pos)+'_'+png_time+'.png'
		png_path = output_dir + path_sep() + png_file
 		write_png,png_path,imb_i,r,g,b
		obs_stc_i.qclk_fname = png_file
		print,'--> Writing 2d observation into a png image file:OK' 
	endif

    if (WRITE_CSV) then begin
        writelog,lunlog,'--> Writing observations output file...'
        obs_file = outfnroot+'_'+cdate_i+'_nan_init.csv'
        obs_path = output_dir + path_sep() + obs_file
        hfc_write_csv,obs_stc_i,obs_path 
        writelog,lunlog,'--> Writing observations output file:OK'
    endif   

    obs_struct = [obs_struct,obs_stc_i]

	writelog,lunlog,'--> Extracting radio sources from the image '+strtrim(i,2)+'...'
	;Extract sources from current image
	feat_stc_i = nrh2d_extract_feat(im_i,ind_i,feat_stc, $
                                    nbmax=in_stc.nbmax, fact=in_stc.fact, $
                                    seuil=in_stc.seuil, fmax=in_stc.fmax, $
                                    amax=in_stc.amax,error=error)
    if (error) then begin
        writelog,lunlog,'WARNING:An error has occured during processing!'
        continue
    endif
 	writelog,lunlog,'--> Extracting radio sources from the image '+strtrim(i,2)+':OK'
   
 	where_feat = where(feat_stc_i.id_rs gt 0L,nfeat)
 	if (where_feat[0] eq -1) then begin
 		writelog,lunlog,'No source found for image '+strtrim(i,2)
 		jd_obs_i = jd_obs_i + jd_cad
		continue
 	endif else writelog,lunlog,strtrim(nfeat,2)+' radio source(s) found.'
    feat_stc_i = feat_stc_i(where_feat)
   
	if (wpng gt 1) then begin
		print,'--> Writing 2d observation + detection results into a png image file...'
		imb_i = bytscl(im_i,top=255,/NAN)
		png_file = file_basename(obs_stc_i.filename)
		png_time = strtrim(strjoin((strsplit(obs_stc_i.date_obs,'T:'))[1:*]),2)
		pos = strpos(png_file,'.',/REVERSE_S)
		png_file = strmid(png_file,0,pos)+'_'+png_time+$
				   '_nrh2d_results.png'
		png_path = output_dir + path_sep() + png_file
		for j=0L,n_elements(feat_stc_i.cc)-1L do begin
			if (strtrim(feat_stc_i(j).cc,2) ne '') or $
			   (strtrim(feat_stc_i(j).cc,2) ne 'NULL') then continue
			cc = feat_cc_extract(feat_stc_i(j).cc,[feat_stc_i(j).cc_x_pix,feat_stc_i(j).cc_y_pix])
			for k=0L,n_elements(cc[0,*])-1L do imb_i[cc[0,k],cc[1,k]] = 255
		endfor
		write_png,png_path,imb_i,r,g,b
		print,'--> Writing 2d observation + detection results into a png image file:OK' 
	endif

 	if (WRITE_CSV) then begin
 		writelog,lunlog,'--> Writing feature parameters output file...'
        feat_file = outfnroot+'_'+cdate_i+'_nan_feat.csv'
        feat_path = output_dir + path_sep() + feat_file
        feat_stc_i.feat_filename = feat_file
        hfc_write_csv,feat_stc_i,feat_path 
  		writelog,lunlog,'--> Writing feature parameters output file:OK'
 	endif
	feat_struct = [feat_struct,feat_stc_i] 

 	if (WRITE_FITS) then begin
		writelog,lunlog,'--> Writing output fits format file...'
        outfname = outfnroot+'_'+cdate_i+'_products.fits'
		nrh2d_write_fits,im_i,index(i),feat_stc_i,$
					 	 write_fits=2,outfname=outfname, $
				     	 output_dir=output_dir,SILENT=SILENT
		writelog,lunlog,'--> Writing output fits format file:OK'
	endif
	jd_obs_i = jd_obs_i + jd_cad
endwhile
writelog,lunlog,'--> Running detection:OK'
if (n_elements(obs_struct) gt 1) then obs_struct = obs_struct(1:*)
if (n_elements(feat_struct) gt 1) then feat_struct = feat_struct(1:*)
;[4]:==============================

writelog,lunlog,'Program has ended correctly on '+systime()
writelog,lunlog,'Elapsed time: '+strtrim(long(systime(1) - syst0),2)+' seconds.'
writelog,lunlog,CLOSE=WRITE_LOG
END