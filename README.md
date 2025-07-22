* Product name: NWLR and PAR
* Responsible person: Hiroshi Murakami (murakami.hiroshi.eo(at)jaxa.jp)
  Research paper: Murakami, Nakayama and Ishizaka 2025, submitted

* Processing time
 Qkm: ~40 min (~26 min with NWLR)
 1km: ~3 min

* Memory
 Qkm: ~2.4 GB
 1km: ~0.2 G

----------------------------

[1] Directory
anc/
bin/
lut/
src/
output/ (not used for the standard operation)

* Temporally files: NA

[2] Input data
 SGLI L1B vnrfile
 SGLI L1B irsfile

[3] Output data
 L2 NWLRfile

[4] Ancillary data (not mandatory):
 ../anc/GGLA/* or ../anc/GGLF/*
 ../anc/ozone/* or ../anc/ozone_jma/*

[5] LUT
 ../lut/
 "values0.dat" and "rgoth_Ebuchi.txt": BRDF correction tables distributed by Morel et al., 2002 through LOV ftp site. (not mandatory)
 LUT_* : atmosphere scattering tables calculated by Pstar4 (Ota et al., 2010)
 vical_pixel.txt: FOV dependency correction table

[6] Return values (PAR)
    0   : success
    1   : not enough execution parameter or Error in h5open_f or h5fclose_f
    2   : No VN-file or No SW-file or No L2-file
    3   : No Scene_center_time or No Latitude or No Longitude or No geometry
    4   : No reflectance
    5   : Error write dsetname
    19  : No parameter
    20  : Compile
    100 : Unexpected error

note:
 * irsfile can be omitted if the resolution/ID is the same as the vnrfile
 * resolution of irsfile can be different from one of vnrfile

---------------------------- Usage -------------------------------
[7] Activation parameters:

 7-0 compile: in O2A_par_pub/bin/
% ./run_pub_sgli_Ariake.sh 0

 7-2 execute: in ./bin/
% ./run_pub_sgli_Ariake.sh $1 ($2) $3 $4

$1: "vnrfile" or 0 (for compile)
$2: "irsfile"
$3: "output-directory" (if it is omitted, it is set to ../output)
$4: AERONET-OC csv file coverd the date of L1B data

---------------------------- Examples ----------------------------
[8] Compile
% cd O2A_par_pub/bin

% ./run_pub_sgli_Ariake i

[9] Execution
% cd O2A_par_pub/bin

test1: output version is to be 2111

% ./run_pub_sgli_Ariake.sh ../input/GC1SG1_202301030209Q06310_1BSG_VNRDQ_2007.h5 ../output ../aeronetoc/ARIAKE_TOWER_L20_2023010120231231.csv
 Aerosol model: AERONET climatology
Open VN file = ../input/GC1SG1_202301030209Q06310_1BSG_VNRDQ_2007.h5
Open SW file = ../input/GC1SG1_202301030209Q06310_1BSG_IRSDQ_2007.h5
 Constant kv =  0.9570 0.9960 0.9760 0.9980 1.0480 1.0340 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
LUCA_error.txt line = 1829
 kt( 1828.1) =  0.9980 1.0009 1.0049 1.0041 1.0045 1.0051 0.9996 1.0017 0.9951 0.9962 0.9970 1.0000 1.0000
 RVS correction: vical_pixel.txt
read: ../aeronetoc/ARIAKE_TOWER_L20_2023010120231231.csv
 -0.0107  33.104 130.272 0.00887 0.01238 0.01486 0.02222 0.02945 0.03545 0.01633 0.01632 0.00581 0.00031 0.00031 0.00031 0.00031 0.00031 0.00031
 Read Land_water_flag
 Read Line_tai93
loc:    5603    5605    3011    3013  -0.011   0.000
 Hokkaido lwfg recovered by MOD44W
 Read TOA refl.: /Image_data/Lt_VN01
 Read TOA refl.: /Image_data/Lt_VN02
 Read TOA refl.: /Image_data/Lt_VN03
 Read TOA refl.: /Image_data/Lt_VN04
 Read TOA refl.: /Image_data/Lt_VN05
 Read TOA refl.: /Image_data/Lt_VN06
 Read TOA refl.: /Image_data/Lt_VN07
 Read TOA refl.: /Image_data/Lt_VN08
 Read TOA refl.: /Image_data/Lt_VN09
 Read TOA refl.: /Image_data/Lt_VN10
 Read TOA refl.: /Image_data/Lt_VN11
 Read TOA refl.: /Image_data/Lt_SW01
 Read TOA refl.: /Image_data/Lt_SW03
read LUT = LUT_pstar4I_SGLI0_lmb_flt_t00_5deg_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o01_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o02_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o03_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o04_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o05_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o06_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o07_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o08_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o09_adj
Read ../anc/GGLA/202301/GGLA_20230103000000.PTW for 2023010300
Read ../anc/GGLA/202301/GGLA_20230103060000.PTW for 2023010306
Read ../anc/GGLA/202301/GGLA_20230103120000.PTW for 2023010312
Read ../anc/GGLA/202301/GGLA_20230103180000.PTW for 2023010318
Read ../anc/GGLA/202301/GGLA_20230104000000.PTW for 2023010400
Read ../anc/ozone/202301/OZONE20230104120000.ctm_le.ana
Read ../anc/ozone/202301/OZONE20230103120000.ctm_le.ana
Read ../anc/ozone/202301/OZONE20230102120000.ctm_le.ana
 Start processing           1
 Proc = 14 / 20   0.17 min
 read = ../lut/rgoth_Ebuchi.txt
 read = ../lut/values0.dat
rsr=  0.1803  0.1799
rs1= 5603 3011  0.0095  0.0133  0.0161  0.0246  0.0330  0.0399  0.0153  0.0153  0.0030  0.0008  0.0008  0.0000  0.0000
rs0= 5603 3011  0.0071  0.0007 -0.0002  0.0115  0.0216  0.0374  0.0153  0.0112 -0.0015  0.0008 -0.0061  0.0040 -0.0010
rar= 5603 3011  0.0146  0.0852  0.1204  0.0984  0.0869  0.0179  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
rsr=  0.1334  0.1330
rs1= 5603 3012  0.0095  0.0133  0.0161  0.0246  0.0330  0.0399  0.0159  0.0159  0.0031  0.0008  0.0008  0.0000  0.0000
rs0= 5603 3012  0.0032  0.0003 -0.0004  0.0127  0.0233  0.0393  0.0159  0.0140  0.0021  0.0008 -0.0039  0.0056  0.0063
rar= 5603 3012  0.0386  0.0865  0.1201  0.0874  0.0720  0.0041  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
rsr=  0.1292  0.1290
rs1= 5604 3012  0.0095  0.0133  0.0161  0.0246  0.0330  0.0399  0.0160  0.0160  0.0031  0.0008  0.0008  0.0000  0.0000
rs0= 5604 3012  0.0097  0.0064  0.0079  0.0157  0.0302  0.0414  0.0160  0.0225  0.0154  0.0008  0.0040  0.0100  0.0038
rar= 5604 3012 -0.0013  0.0484  0.0622  0.0699  0.0214 -0.0110  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000
 Start processing           2
 Proc =  0 / 20   0.17 min
 Proc =  1 / 20   0.48 min
 Proc =  2 / 20   0.78 min
 Proc =  3 / 20   1.21 min
 Proc =  4 / 20   1.85 min
 Proc =  5 / 20   2.57 min
 Proc =  6 / 20   3.29 min
 Proc =  7 / 20   4.06 min
 Proc =  8 / 20   4.93 min
 Proc =  9 / 20   5.77 min
 Proc = 10 / 20   6.54 min
 Proc = 11 / 20   7.14 min
 Proc = 12 / 20   7.53 min
 Proc = 13 / 20   7.97 min
 Proc = 14 / 20   8.53 min
rc1= 5604 3012  0.0738  0.0820  0.0849  0.0892  0.0883  0.0855  0.0719  0.0768  0.0413  0.0590  0.0616  0.0624  0.0498
rts= 5604 3012  0.0732  0.0808  0.0834  0.0862  0.0852  0.0828  0.0719  0.0719  0.0368  0.0590  0.0590  0.0543  0.0469
  5604  3012  0.0124  0.0124  1.0734
 Proc = 15 / 20   9.10 min
 Proc = 16 / 20   9.69 min
 Proc = 17 / 20  10.31 min
 Proc = 18 / 20  10.93 min
 Proc = 19 / 20  11.55 min
 Proc = 20 / 20  12.14 min
Create L2 file = ../output/GC1SG1_202301030209Q06310_L2SG_NWLRQ_2111.h5
 r_out = 0.65752
 Get processing date
Write /Image_data/PAR
Write /Image_data/NWLR_380
Write /Image_data/NWLR_412
Write /Image_data/NWLR_443
Write /Image_data/NWLR_490
Write /Image_data/NWLR_530
Write /Image_data/NWLR_565
Write /Image_data/NWLR_670
Write /Image_data/CHLA
Write /Image_data/aph_442
Write /Image_data/CDOM
Write /Image_data/TSM
Write /Image_data/TAUA_865
Write /Image_data/TAUA_670
Write /Image_data/FAI
Write /Image_data/QA_flag
Write /Image_data/Line_tai93
 Deallocate Arrays
 Finish
Proc time: 20250718 00:09:08 - 20250718 00:21:40
ret = 0

---------
test2 (aeronet OC is not available on the data day): output version is to be 2011

./run_pub_sgli_Ariake.sh /export/bass/archive/eorc/GCOM-C/original/SGLI/L1SG/2/1B_/201907/26/GC1SG1_201907260144B05410_1BSG_VNRDQ_2000.h5 ../outpWER_L20_2019010120191231.csv
AlmaLinux release 8.4 (Electric Cheetah) on kona81, exe: pub_sgli_Ariake_8i
% run_pub_sgli_Ariake.sh /export/bass/archive/eorc/GCOM-C/original/SGLI/L1SG/2/1B_/201907/26/GC1SG1_201907260144B05410_1BSG_VNRDQ_2000.h5 ../output ../aerone010120191231.csv
 Aerosol model: AERONET climatology
Open VN file = ../input/GC1SG1_201907260144B05410_1BSG_VNRDQ_2000.h5
Open SW file = ../input/GC1SG1_201907260144B05410_1BSG_IRSDQ_2000.h5
 Constant kv =  0.9570 0.9960 0.9760 0.9980 1.0480 1.0340 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
LUCA_error.txt line =  572
 kt(  571.1) =  1.0004 0.9997 1.0007 1.0002 1.0000 0.9996 0.9994 0.9990 0.9998 0.9996 0.9999 1.0000 1.0000
 RVS correction: vical_pixel.txt
read: ../aeronetoc/ARIAKE_TOWER_L20_2019010120191231.csv
 Read Land_water_flag
 Read Line_tai93
loc:    5989    5991     817     819   0.000   0.001
 Hokkaido lwfg recovered by MOD44W
 Read TOA refl.: /Image_data/Lt_VN01
 Read TOA refl.: /Image_data/Lt_VN02
 Read TOA refl.: /Image_data/Lt_VN03
 Read TOA refl.: /Image_data/Lt_VN04
 Read TOA refl.: /Image_data/Lt_VN05
 Read TOA refl.: /Image_data/Lt_VN06
 Read TOA refl.: /Image_data/Lt_VN07
 Read TOA refl.: /Image_data/Lt_VN08
 Read TOA refl.: /Image_data/Lt_VN09
 Read TOA refl.: /Image_data/Lt_VN10
 Read TOA refl.: /Image_data/Lt_VN11
 Read TOA refl.: /Image_data/Lt_SW01
 Read TOA refl.: /Image_data/Lt_SW03
read LUT = LUT_pstar4I_SGLI0_lmb_flt_t00_5deg_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o01_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o02_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o03_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o04_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o05_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o06_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o07_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o08_adj
read LUT = LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o09_adj
No ANC data : ../anc/GGLF/201907/GGLF_20190726000000.PTW
No ANC data : ../anc/GGLF/201907/GGLF_20190726060000.PTW
No ANC data : ../anc/GGLF/201907/GGLF_20190726120000.PTW
No ANC data : ../anc/GGLF/201907/GGLF_20190726180000.PTW
No ANC data : ../anc/GGLF/201907/GGLF_20190727000000.PTW
 No ozone data: = 343.79 DU, day:  1
 No ozone data: = 343.79 DU, day:  2
 No ozone data: = 343.79 DU, day:  3
 Start processing           1
 Proc =  0 / 20   0.16 min
 No file: lut/rgoth_Ebuchi.txt
 Cannot read: ../lut/values0.dat
 Proc =  1 / 20   0.50 min
 Proc =  2 / 20   0.85 min
 Proc =  3 / 20   1.21 min
 Proc =  4 / 20   1.70 min
 Proc =  5 / 20   2.27 min
 Proc =  6 / 20   2.79 min
 Proc =  7 / 20   3.36 min
 Proc =  8 / 20   3.98 min
 Proc =  9 / 20   4.60 min
 Proc = 10 / 20   5.18 min
 Proc = 11 / 20   5.65 min
 Proc = 12 / 20   6.10 min
 Proc = 13 / 20   6.50 min
 Proc = 14 / 20   7.00 min
 Proc = 15 / 20   7.59 min
 Proc = 16 / 20   8.18 min
 Proc = 17 / 20   8.80 min
 Proc = 18 / 20   9.48 min
 Proc = 19 / 20  10.19 min
 Proc = 20 / 20  10.90 min
Create L2 file = ../output/GC1SG1_201907260144B05410_L2SG_NWLRQ_2011.h5
 r_out = 0.75191
 Get processing date
Write /Image_data/PAR
Write /Image_data/NWLR_380
Write /Image_data/NWLR_412
Write /Image_data/NWLR_443
Write /Image_data/NWLR_490
Write /Image_data/NWLR_530
Write /Image_data/NWLR_565
Write /Image_data/NWLR_670
Write /Image_data/CHLA
Write /Image_data/aph_442
Write /Image_data/CDOM
Write /Image_data/TSM
Write /Image_data/TAUA_865
Write /Image_data/TAUA_670
Write /Image_data/FAI
Write /Image_data/QA_flag
Write /Image_data/Line_tai93
 Deallocate Arrays
 Finish
Proc time: 20250718 00:02:15 - 20250718 00:13:31
ret = 0
