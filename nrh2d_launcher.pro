PRO nrh2d_launcher,starttime,endtime,$
				    cadence=cadence,hend=hend,$
				    config_file=config_file,$
					freq=freq,data_dir=data_dir,$
					script_dir=script_dir,$
					output_dir=output_dir,$
					config_dir=config_dir,$
					write_png=write_png, $ 
					WRITE_CSV=WRITE_CSV,$
					WRITE_FITS=WRITE_FITS,$
					CLEAN_DATA=CLEAN_DATA,$
					DOWNLOAD_DATA=DOWNLOAD_DATA,$
                    GET_CLOSEST_FREQ=GET_CLOSEST_FREQ, $
					SILENT=SILENT

;+
; NAME:
;		nrh2d
;
; PURPOSE:
; 		Launcher for the nrh2d software. 
;
; CATEGORY:
;		Image processing
;
; GROUP:
;		NRH2D
;
; CALLING SEQUENCE:
;		IDL>nrh2d_launcher,starttime,endtime
;
; INPUTS:
;		starttime   - First date of time range to process.
;					  (Format is YYYY-MM-DD.)
;		endtime     - Last date of time range to process.
;					  (Format is YYYY-MM-DD.) 
;	
; OPTIONAL INPUTS:
;		cadence     - Cadence of the images in seconds.
;					  Default is 10 seconds.
;		hend	    - Scalar of string type containing the 
;				 	  time of the last image to process for
;					  each day (in UTC).
;					  Default is '15:00'
;
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
;					      Defaults is 150 MHz.
;		config_file - Scalar of string type containing the name of the configuration file.
;		config_dir  - Scalar of string type providing the path to the configuration file directory. 
;		data_dir    - Scalar of string type providing the path to the input nrh data file directory.
;		script_dir  - Scalar of string type providing the path to the scripts directory.
;		output_dir  - Scalar of string type containing the path to the 
;					      directory where output data files will be saved.
;					      (If output_dir is not specified, then the program
;					      will create a directory Products in the current 
;					      folder.)
;		write_png   - Write output png image files:
;							write_png = 0 --> No output
;							write_png = 1 --> Write png 2d image of the radio observation
;											  without detection results.
;							write_png = 2 --> Write png 2d imageof the radio obervation
;											  including detection results.
;
; KEYWORD PARAMETERS:
;		/DOWNLOAD_DATA    - Allow program to download data
;						    from distant ftp server (if available).
;		/CLEAN_DATA	      - Remove data file after process.
;		/WRITE_CVS        - Write cvs files containing detection products.
;		/WRITE_FITS	      - Write fits file of the 2d image with detection
;					        results overplotted.
;       /GET_CLOSEST_FREQ - If the input frequency is not available, get the closet one available.
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
; RESTRICTIONS/COMMENTS:
;		The Solar SoftWare (SSW) with the NRH library package 
;		must be loaded.
;			
; CALL:
;		countday
;		nrh2d
;       getftp_nrh10s.csh (script)
;
; EXAMPLE:
;		None.		
;
; MODIFICATION HISTORY:
;		Written by X.Bonnin,	26-SEP-2011.			
;				
;-

if (n_params() lt 2) then begin
	message,/INFO,'Call is:'
	print,'nrh2d_launcher,starttime,endtime,$'
	print,'                config_file=config_file,$'
	print,'                freq=freq,data_dir=data_dir,$'
	print,'                cadence=cadence,hend=hend,$'
	print,'                output_dir=output_dir,$'
	print,'                script_dir=script_dir,config_dir=config_dir,$'
	print,'                write_png=write_png,/DOWNLOAD_DATA,/CLEAN_DATA,$'
	print,'                /SILENT,/WRITE_CSV,/WRITE_FITS,/GET_CLOSEST_FREQ'
	return
endif

DOWNLOAD_DATA = keyword_set(DOWNLOAD_DATA)
SILENT = keyword_set(SILENT)
WRITE_CSV = keyword_set(WRITE_CSV)
WRITE_FITS=keyword_set(WRITE_FITS)
GET_CLOSEST_FREQ = keyword_set(GET_CLOSEST_FREQ)
CLEAN_DATA = keyword_set(CLEAN_DATA)

startt = strjoin(strsplit(strtrim(starttime[0],2),'-',/EXTRACT))
endt = strjoin(strsplit(strtrim(endtime[0],2),'-',/EXTRACT))
date = countday(startt,endt,nday=nday)

cd,current=curdir
if (~keyword_set(data_dir)) then datdir=curdir else datdir = strtrim(data_dir[0],2)
if (~keyword_set(output_dir)) then outdir=curdir else outdir = strtrim(output_dir[0],2)
if (~keyword_set(config_file)) then begin
    ifile=file_search('../','nrh2d_*_config.txt',count=nifile) 
    if (nifile ne 1) then begin
        message,/CONT,'configuration file not found!'
        return
    endif
    ifile = ifile[0]
endif else ifile = strtrim(config_file[0],2)
if (~keyword_set(config_dir)) then indir = file_dirname(ifile) else indir = strtrim(config_dir[0],2)
if (~keyword_set(script_dir)) then scptdir = curdir else scptdir = strtrim(script_dir[0],2) 

file = '2i'+strmid(date,2)+'.*'

for i=0L,nday-1L do begin

	if (DOWNLOAD_DATA) then begin
		yyyy = strmid(date[i],0,4)
		mm = strmid(date[i],4,2)
		dd = strmid(date[i],6,2)
		file_i = file_search(datdir + path_sep() + file[i],count=nfile_i)
		if (file_i[0] eq '') then begin
            if not (SILENT) then print,'downloading '+file[i]+'...'
			cd,datdir
			spawn,'csh '+scptdir+path_sep()+'getftp_nrh10s.csh '+dd+' '+mm+' '+yyyy
			cd,curdir
		endif
	endif

	file_i = file_search(datdir + path_sep() + file[i],count=nfile_i)
	if (file_i[0] eq '') then continue

	for j=0l, nfile_i-1l do begin
       if not (SILENT) then print,'Processing '+file_i[j]+'...'

	   nrh2d,file_i[j],ifile,$
	  	    freq=freq,input_dir=indir,$
	  	    cadence=cadence,hend=hend,$
	   	    output_dir=outdir,$
	 	    write_png=write_png, $
	   	    WRITE_FITS=WRITE_FITS,$
	   	    WRITE_CSV=WRITE_CSV,$
	   	    GET_CLOSEST_FREQ=GET_CLOSEST_FREQ, $
            SILENT=SILENT
       if not (SILENT) then print,'Processing '+file_i[j]+':OK' 

	   if (CLEAN_DATA) then begin
            if not (SILENT) then print,'Deleting '+file_i[j]+'...'
            spawn,'rm -f '+file_i[j]
       endif
    endfor
endfor

END