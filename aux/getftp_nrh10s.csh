#! /bin/csh

if ($# == 3) then
    set month = ${2}
    set day = ${1}
    set year_yy = `echo ${3} | cut -c3-4`
    set year_yyyy = ${3}
else
    echo "Call is: csh getftp_nrh10s day month year"
    exit
endif

if ($year_yyyy > 1998) then
	set data_dir = "/data2/juke/10sec"
else 
	set data_dir = "/data2/juke/10secbis"
endif

if ($day <= 10) set rep = "A"
if ($day > 10 && $day <= 20) set rep = "B"	
if ($day > 20 && $day <= 31) set rep = "C"

# FTP pour reccuperer les fichiers NRH 10 sec
# ATTENTION au binary avec mesopl nouveau sinon les fichiers sont corrompus
ftp -nv <<%
open mesopl.obspm.fr
user bonnin Mesopl@260183!
cd $data_dir"/INT"$year_yy$month$rep"/"
binary
mget "2i"$year_yy$month$day".*"
a
close
%

# Sert a de pas avoir d'erreur du style No match
set nonomatch

exit
