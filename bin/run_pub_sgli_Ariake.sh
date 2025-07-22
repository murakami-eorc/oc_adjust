#!/bin/csh -f
#
set hname0 = `hostname`
set osvsn0 = `cat /etc/redhat-release`
#
set progfile="../src/pub_sgli_v209_Ariake.f"
#
#------------------------------------- HDF5-env
set cpl = `ls -t1 ./pub_sgli_Ariake_8* | head -1 | tail -c2`
if ( $1 == "i" || $1 == "-i" || ( $1 != "0" && $1 != "-1" && $cpl == "i" ) ) then
 set binfile="pub_sgli_Ariake_8i"
 set PrefixPATH=/export/emc2/util_RHEL8/intel/hdf5-1.8.21_disable_szip
 set SZPATH=/export/emc2/util_RHEL8/intel/szip-2.1.1
 set ZPATH=/export/emc2/util_RHEL8/intel/zlib-1.2.8
else
 set binfile="pub_sgli_Ariake_8p"
 set PrefixPATH=/export/emc2/util_RHEL8/pgi/hdf5-1.8.21_disable_szip
 set SZPATH=/export/emc2/util_RHEL8/pgi/szip-2.1.1
 set ZPATH=/export/emc2/util_RHEL8/gnu/zlib-1.2.11
endif
echo $osvsn0" on "$hname0", exe: "$binfile
#
if ( ( $1 == "0" ) || ( $1 == "-1" ) || ( $1 == "i" ) || ( $1 == "-i" ) || ( $1 == "g" ) ) then
 if (! $?LD_LIBRARY_PATH) then 
  setenv LD_LIBRARY_PATH $PrefixPATH/lib
 else
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$PrefixPATH/lib
 endif
 if (! $?LD_RUN_PATH) then 
  setenv LD_RUN_PATH $PrefixPATH/lib
 else
  setenv LD_RUN_PATH ${LD_RUN_PATH}:$PrefixPATH/lib
 endif
 setenv HDF_INC $PrefixPATH/include
 setenv HDF_LIB $PrefixPATH/lib
 setenv PATH    $PrefixPATH/bin:$PATH
endif
#---------------------------------------------
#
if (-e "../bin/")  then
 set bindir="../bin/"
else
 set bindir="./"
endif
#
if ($1 == "") then
 echo "> run_pub_sgli_Ariake.sh vnfile swfile outdir(or outfile) (AERONET-OC.csv)"
 exit 19
else if ($1 == "0") then
 pgf95 -O3 -Bstatic $progfile -o $bindir$binfile \
       -I$HDF_INC -L$HDF_LIB -lhdf5_fortran -lhdf5 -lz\
       -I$SZPATH/include -L$SZPATH/lib \
       -I$ZPATH/include -L$ZPATH/lib
 exit 20
else if ($1 == "-1") then
 pgf95 -g -Ktrap=fp -C -Bstatic $progfile -o $bindir$binfile \
       -I$HDF_INC -L$HDF_LIB -lhdf5_fortran -lhdf5 -lz\
       -I$SZPATH/include -L$SZPATH/lib \
       -I$ZPATH/include -L$ZPATH/lib
 exit 20
else if ($1 == "i") then
 ifort -O3 -Bstatic -assume bscc -assume byterecl -ftz -ip -ipo $progfile -o $bindir$binfile\
       -I$HDF_INC -L$HDF_LIB -lhdf5_fortran -lhdf5 -lz \
       -I$SZPATH/include -L$SZPATH/lib\
       -I$ZPATH/include -L$ZPATH/lib
 exit 20
else if ($1 == "-i") then
 ifort -g -assume byterecl -check all -warn all $progfile -o $bindir$binfile\
       -I$HDF_INC -L$HDF_LIB -lhdf5_fortran -lhdf5 -lz \
       -I$SZPATH/include -L$SZPATH/lib\
       -I$ZPATH/include -L$ZPATH/lib
 exit 20
else
 echo "% run_pub_sgli_Ariake.sh "$1 $2 $3 $4
 set dates=`date -u '+%Y%m%d %T'`
# pgdbg ./$binfile $1 $2 $3 $4
 ./$binfile $1 $2 $3 $4
 set ret0 = $status
 set datee=`date -u '+%Y%m%d %T'`
 echo "Proc time: "$dates' - '$datee
#
 if ($ret0 == 100) then
  set ret = 0
 else if ($ret0 == 0) then
  set ret = 100
 else
  set ret = $ret0
 endif
 echo "ret = "$ret
 exit $ret
endif
