echo -e "\
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   gmx_mmpbsa   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    Jicun Li    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
>>>>>>>>>>>>>>>>>>>>>>>>>>     2022-11-18 00:24:41     <<<<<<<<<<<<<<<<<<<<<<<<<\n
>>   Usage: gmx_mmpbsa [-f   *.xtc  ] [-s   *.tpr  ] [-n   *.ndx     ]
                       [-pro PROTEIN] [-lig <LIGAND|none>]
                       [-cou dh     ] [-ts ie      ] [-cas <B:E:S,#> ]\n
>> Default: gmx_mmpbsa -f   traj.xtc  -s   topol.tpr  -n   index.ndx
                       -pro Protein   -lig Ligand
                                       -ts ie
--------------------------------------------------------------------------------
>> Option:
      -f: trajectory file
      -s: topology file
      -n: index file
    -pro: index group name of protein
    -lig: index group name of ligand, can be ignored using none
    -cou: dh: calculate MM(COU) with Debye-Huckel screening method
     -ts: ie: calculate entropy with interaction entropy method
    -cas: select res for CAS, begin:end:step or number, using global resnum
   Other: change them in the script directly
--------------------------------------------------------------------------------
>> Warning:
   !!! 1. make sure gmx is available and version is consistent
   !!! 2. gawk version > 3; macOS awk version > 2018
   !!! 3. watch _pid.pdb to make sure PBC is correct
--------------------------------------------------------------------------------
>> Log:
   TODO:       parallel APBS, focus
   2022-11-18: fix issues with non-continuous ndx
   2022-10-16: fix issue with ndx file
   2022-07-14: fix issues with Chinease in PBSAset
   2022-07-08: fix bug, No PB SA value
   2022-07-01: fix bug for awk too many open files
   2022-06-18: fix bug for ndx when Lig is 1st
   2022-04-19: fix bug for CAS
               remove -com option
   2022-02-09: CAS except GLY, PRO, CYX
               fix radius bug, use AMBER mBondi
               change to AMBER PB4
   2021-11-25: keep original PDB residue name
               output averaged contributions
               clean output files
   2021-05-17: add systime and fflush, only for gawk
   2021-05-16: revise cmd for all Linux
               see https://github.com/Jerkwin/gmxtool/issues/2
   2021-03-14: add -cou, -ts option, see 10.1088/0256-307X/38/1/018701
               change default PB method to npbe
   2021-03-03: fix bug: awk will exit for multiple system call
   2020-11-15: fix bug for resID >=1000, for awk 3.x gsub /\s/
   2020-06-02: fix bug of withLig
   2020-06-01: fix bug of -skip
   2020-05-27: fix bug for sdie
   2020-05-26: fix bug for RES name
   2020-04-03: use C6, C12 directly
   2020-01-08: support protein only
   2019-12-24: fix bug for small time step
   2019-12-10: fix bug for OPLS force field
   2019-11-17: fix bug for c6, c12 of old version tpr
   2019-11-03: fix bug for time stamp
   2019-09-19: push to gmxtool
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"

################################################################################
# 设置运行环境, 计算参数
# setting up environments and parameters
################################################################################

trj=traj.xtc		# 轨迹文件 trajectory file
tpr=topol.tpr		# tpr文件  tpr file
ndx=index.ndx		# 索引文件 index file

pro=Protein			# 蛋白索引组   index group name of protein
lig=Ligand			# 配体索引组   index group name of ligand

step=0	# 从第几步开始运行 step number to run
		# 1. 预处理轨迹: 复合物完整化, 团簇化, 居中叠合, 然后生成pdb文件
		# 2. 获取每个原子的电荷, 半径, LJ参数, 然后生成qrv文件
		# 3. MM-PBSA计算: pdb->pqr, 输出apbs, 计算MM, APBS
		# 1. pre-processe trajectory, whole, cluster, center, fit, then generate pdb file
		# 2. abstract atomic parameters, charge, radius, C6/C12, then generate qrv file
		# 3. run MM-PBSA, pdb->pqr, apbs, then calculate MM, PB, SA

gmx='gmx'								# /path/to/GMX/bin/gmx_mpi
dump="$gmx dump"						# gmx dump
trjconv="$gmx trjconv -dt 1000"			# gmx trjconv, use -b -e -dt, NOT -skip
converttpr="$gmx convert-tpr"			# to generate reduced tpr

#apbs='/path/APBS/bin/apbs'				# APBS(Linux)
apbs='C:/Software/APBS/bin/apbs.exe'				# APBS(Windows), USE "/", NOT "\"
export MCSH_HOME=/dev/null				# APBS io.mc

pid=pid				# 输出文件($$可免重复) prefix of the output files($$)
err=_$pid.err		# 屏幕输出文件 file to save the message from the screen
qrv=_$pid.qrv		# 电荷/半径/VDW参数文件 to save charge/radius/vdw parameters
pdb=_$pid.pdb		# 轨迹构型文件 to save trajectory
idx=_$pid.ndx		# 临时索引文件 to generate pdb
tpx=_$pid.tpr		# 临时tpr文件  to generate pdb

radType=1			# 原子半径类型 radius of atoms (0:ff; 1:mBondi)
radLJ0=1.2			# 用于LJ参数原子的默认半径(A, 主要为H) radius when LJ=0 (H)

meshType=0			# 网格大小 mesh (0:global  1:local)
gridType=1			# 格点大小 grid (0:GMXPBSA 1:psize)

cfac=3				# 分子尺寸到粗略格点的放大系数
					# Factor to expand mol-dim to get coarse grid dim
fadd=10				# 分子尺寸到细密格点的增加值(A)
					# Amount added to mol-dim to get fine grid dim (A)
df=0.5				# 细密格点间距(A) The desired fine mesh spacing (A)

# 极性计算设置(Polar)
PBEset='
  temp  298.15      # 温度
  pdie  2           # 溶质介电常数
  sdie  78.54       # 溶剂介电常数, 真空1, 水78.54

  npbe              # PB方程求解方法, lpbe(线性), npbe(非线性), smbpe(大小修正)
  bcfl  mdh         # 粗略格点PB方程的边界条件, zero, sdh/mdh(single/multiple Debye-Huckel), focus, map
  srfm  smol        # 构建介质和离子边界的模型, mol(分子表面), smol(平滑分子表面), spl2/4(三次样条/7阶多项式)
  chgm  spl4        # 电荷映射到格点的方法, spl0/2/4, 三线性插值, 立方/四次B样条离散
  swin  0.3         # 立方样条的窗口值, 仅用于 srfm=spl2/4

  srad  1.4         # 溶剂探测半径
  sdens 10          # 表面密度, 每A^2的格点数, (srad=0)或(srfm=spl2/4)时不使用

  ion charge  1 conc 0.15 radius 0.95  # 阳离子的电荷, 浓度, 半径
  ion charge -1 conc 0.15 radius 1.81  # 阴离子

  calcforce  no
  calcenergy comps'

# 非极性计算设置(Apolar/Non-polar)
PBAset='
  temp  298.15 # 温度
  srfm  sacc   # 构建溶剂相关表面或体积的模型
  swin  0.3    # 立方样条窗口(A), 用于定义样条表面

  # SASA
  srad  1.4    # 探测半径(A)
  gamma 1      # 表面张力(kJ/mol-A^2)

  #gamma const 0.0226778 3.84928
  #gamma const 0.027     0
  #gamma const 0.0301248 0         # AMBER-PB4 .0072*cal2J 表面张力, 常数

  press  0     # 压力(kJ/mol-A^3)
  bconc  0     # 溶剂本体密度(A^3)
  sdens 10
  dpos  0.2
  grid  0.1 0.1 0.1

  # SAV
  #srad  1.29      # SAV探测半径(A)
  #press 0.234304  # 压力(kJ/mol-A^3)

  # WCA
  #srad   1.25           # 探测半径(A)
  #sdens  200            # 表面的格点密度(1/A)
  #dpos   0.05           # 表面积导数的计算步长
  #bconc  0.033428       # 溶剂本体密度(A^3)
  #grid   0.45 0.45 0.45 # 算体积分时的格点间距(A)

  calcforce no
  calcenergy total'

################################################################################
# 检查 gmx, apbs 是否可以运行
# check gmx, apbs availability
################################################################################

str=$($gmx --version | grep -i "GROMACS version")
[[ $step -lt 3 && -z "$str" ]] && { echo -e "!!! ERROR !!!  GROMACS NOT available !\n"; exit; }

str=$($apbs --version 2>&1 | grep -i "Poisson-Boltzmann")
[[ -z "$str" ]] && { echo -e "!!! WARNING !!!  APBS NOT available !\n"; }

################################################################################
# 解析命令行参数
# parse command line options
################################################################################

useDH=1;  useTS=1; isCAS=0

cas="$*";
cas=${cas##*-cas}; cas=${cas%%-*}; cas=${cas// /}

opt=($*); N=${#opt[@]}
for((i=0; i<N; i++)); do
	arg=${opt[$i]}; j=$((i+1)); val=${opt[$j]}
	[[ $arg =~ -f   ]] && { trj=$val; }
	[[ $arg =~ -s   ]] && { tpr=$val; }
	[[ $arg =~ -n   ]] && { ndx=$val; }
	[[ $arg =~ -pro ]] && { pro=$val; }
	[[ $arg =~ -lig ]] && { lig=$val; }
	[[ $arg =~ -cou ]] && { useDH=1;  }
	[[ $arg =~ -ts  ]] && { useTS=1;  }
	[[ $arg =~ -cas ]] && { isCAS=1;  }
done

if [[ -z "$lig" || $lig =~ none ]]; then
	withLig=0; com=$pro; lig=$pro;
else
	withLig=1; com=${pro}_${lig}
fi

if [[ $step -lt 3 ]]; then
################################################################################
# 检查所需文件, 索引组
# check needed files, index group
################################################################################

[[ ! -f "$trj" ]] && { echo -e "!!! ERROR !!! trajectory File NOT Exist !\n"; exit; }
[[ ! -f "$tpr" ]] && { echo -e "!!! ERROR !!! topology   File NOT Exist !\n"; exit; }
[[ ! -f "$ndx" ]] && { echo -e "!!! ERROR !!! index      File NOT Exist !\n"; exit; }

str=$(grep "\[ *$pro *\]" "$ndx"); [[ -z "$str"                   ]] && { echo -e "!!! ERROR !!! [ $pro ] NOT in $ndx !\n"; exit; }
str=$(grep "\[ *$lig *\]" "$ndx"); [[ -z "$str" && $withLig -eq 1 ]] && { echo -e "!!! ERROR !!! [ $lig ] NOT in $ndx !\n"; exit; }
str=$(grep "\[ *$pro *\]" "$ndx" | wc -l); [[ $str -gt 1                   ]] && { echo -e "!!! ERROR !!! more than ONE [ $pro ] in $ndx !\n"; exit; }
str=$(grep "\[ *$lig *\]" "$ndx" | wc -l); [[ $str -gt 1 && $withLig -eq 1 ]] && { echo -e "!!! ERROR !!! more than ONE [ $lig ] in $ndx !\n"; exit; }

if [[ $withLig -eq 1 ]]; then
	awk -v pro=$pro -v lig=$lig '
		BEGIN {i=-1; cp=""}
		index($0, "[") { i++
			txt=$0;
			sub(/ *\[ */, "", txt);
			sub(/ *\].*/, "", txt)
			if(txt==pro || txt==lig) cp=cp" "i
		}
		END {
			sub(/^ */, "", cp)
			gsub(/ +/, "\n", cp)
			print cp"\ndel 0-"i"\n0|1\nname 2 "pro"_"lig"\nq"
		}
	' $ndx | $gmx make_ndx -f $tpr -n $ndx -o $idx &> $err
fi

echo -e ">> 0. set up environments and parameters: OK !\n"
fi

if [[ $step -le 1 ]]; then
################################################################################
# 1. 预处理轨迹: 复合物完整化, 团簇化, 居中叠合, 然后生成pdb文件
#    请检查pdb文件确保构型PBC处理正确
# 1. pre-processe trajectory, whole, cluster, center, fit, then generate pdb file
################################################################################

trjwho=_$pid~trj_who.xtc;
echo "$com" | $trjconv -s $tpr -n $idx -f $trj -o $trjwho -pbc whole &>>$err
[[ ! -f "$trjwho" ]] && { echo -e "!!! ERROR !!! gmx trjconv Failed ! Check $err to figure out why.\n"; exit; }

echo "$com" | $converttpr -s $tpr -n $idx -o $tpx &>>$err

awk >$idx -v ndx=$ndx -v pro=$pro -v lig=$lig -v withLig=$withLig '
	BEGIN { RS="["
		while(getline < ndx) {
			gsub(" ","", $1); gsub("\t","", $1)
			if($1==pro)    { ipro=$3; npro=NF-2; }
			if($1==pro"]") { ipro=$2; npro=NF-1; }
			if(withLig) {
				if($1==lig)    { ilig=$3; nlig=NF-2; }
				if($1==lig"]") { ilig=$2; nlig=NF-1; }
			}
		}
		if(ipro<ilig) { dpro=0; dlig=npro }
		else          { dlig=0; dpro=nlig }
		print "[ "pro" ]";       for(i=1; i<=npro; i++)      printf "%d%s", i+dpro, i%15==0?"\n":" "; print ""
		print "[ "lig" ]";       for(i=1; i<=nlig; i++)      printf "%d%s", i+dlig, i%15==0?"\n":" "; print ""
		print "[ "pro"_"lig" ]"; for(i=1; i<=npro+nlig; i++) printf "%d%s", i,      i%15==0?"\n":" "; print ""
	} '

trjcnt=_$pid~trj_cnt.xtc; trjcls=_$pid~trj_cls.xtc
if [[ $withLig -eq 1 ]]; then
# usful for single protein and ligand
	echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $trjcnt &>>$err -pbc mol -center
	echo -e "$com\n$com" | $trjconv  -s $tpx -n $idx -f $trjcnt -o $trjcls &>>$err -pbc cluster
	echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjcls -o $pdb    &>>$err -fit rot+trans
else
	echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $pdb    &>>$err -pbc mol -center
fi
echo -e ">> 1. preprocessing trajectory: OK !\n"

fi; if [[ $step -le 2 ]]; then
################################################################################
# 2. 获取每个原子的电荷, 半径, LJ参数, 然后生成qrv文件
# 2. abstract atomic parameters, charge, radius, C6/C12, then generate qrv file
#    feel free to change radius with radType
#    radType=0: radius from C6/C12, or radLJ0 if either C6 or C12 is zero
#    radType=1: mBondi
################################################################################

$dump -quiet -s $tpr 2>>$err \
| awk >$qrv -v ndx=$ndx -v pro=$pro -v lig=$lig -v withLig=$withLig \
			-v radType=$radType -v radLJ0=$radLJ0 '
	BEGIN { RS="["
		print pro, lig
		while(getline < ndx) {
			gsub(" ","", $1); gsub("\t","", $1)
			if($1==pro)    for(i=3; i<=NF; i++) ndxPro[$i+0]++
			if($1==pro"]") for(i=2; i<=NF; i++) ndxPro[$i+0]++
			if(withLig) {
				if($1==lig)    for(i=3; i<=NF; i++) ndxLig[$i+0]++
				if($1==lig"]") for(i=2; i<=NF; i++) ndxLig[$i+0]++
			}
		}
		RS="\r?\n"
		nres=0
	}

	/#molblock/  { Ntyp=$3 }
	/moltype.+=/ { Imol=$3; getline; Nmol[Imol]=$3 }
	/ffparams:/ {
		getline Atyp; sub(/.+=/, "", Atyp); Atyp += 0
		print Atyp
		getline
		for(i=0; i<Atyp; i++) {
			printf "%6d", i
			for(j=0; j<Atyp; j++) {
				getline
				C6 =$0; sub(".*c6 *= *",  "", C6);  sub(",.*", "", C6);
				C12=$0; sub(".*c12 *= *", "", C12); sub(",.*", "", C12);
				printf " %s %s", C6, C12
				if(j==i) {
					sigma[i]=0; epsilon[i]=0
					Rad[i]=radLJ0
					if(C6*C12!=0) {
						sigma[i]=10*(C12/C6)^(1./6) # 转换单位为A
						epsilon[i]=C6^2/(4*C12)
						Rad[i]=.5*sigma[i]          # sigma为直径
					}
				}
			}
			print ""
		}
	}

	/moltype.+\(/ { Imol=$0; gsub(/[^0-9]/,"",Imol)
		getline txt; sub(/.*=/,"",txt); gsub(" ","_",txt)
		Name[Imol]=txt
		getline; getline txt;       gsub(/[^0-9]/,"",txt); Natm[Imol]=txt+0
		for(i=0; i<Natm[Imol]; i++) {
			getline; txt=$0; idx=$3; resID[Imol, i]=$(NF-2)+1+nres
			sub(",", "", idx);    idx += 0;
			Catm[Imol, i]=idx
			Ratm[Imol, i]=Rad[idx]
			Satm[Imol, i]=sigma[idx]
			Eatm[Imol, i]=epsilon[idx]
			sub(/.+q=/, "", txt); sub(/,.+/,  "", txt); Qatm[Imol, i]=txt
		}
		getline
		for(i=0; i<Natm[Imol]; i++) {
			getline txt
			sub(/.+=./, "", txt); sub(/..$/, "", txt)
			Tatm[Imol, i]=txt
		}
	}

	/residue\[/ { nres++
		txt=$0; sub(/.*nr=/,"",txt); sub(/,.*/,"",txt); txt*=1
		sub(/.*="/,"",$0); sub(/".*/,"",$0);
		resName[nres]=sprintf("%05d%s", txt, $0)
	}

	/^ +Angle:/ {
		getline; n=$2+0
		if(n>0) {
			getline
			while(getline) {
				if(!index($0, "ANGLES")) break
				k=$NF; j=$(NF-1); i=$(NF-2)
				tag=toupper(Tatm[Imol,i])
				if(substr(tag,1,1)=="H") Tatm[Imol,i]=sprintf("H%s", Tatm[Imol,j])
				tag=toupper(Tatm[Imol,k])
				if(substr(tag,1,1)=="H") Tatm[Imol,k]=sprintf("H%s", Tatm[Imol,j])
			}
		}
	}

	END {
		Ntot=0; Nidx=0
		for(i=0; i<Ntyp; i++) {
			for(n=0; n<Nmol[i]; n++) {
				for(j=0; j<Natm[i]; j++) {
					Ntot++
					if(Ntot in ndxPro || Ntot in ndxLig) {
						Nidx++
						if(radType==0) radi=Ratm[i, j]
						if(radType >0) radi=getRadi(Tatm[i, j], radType)
						printf "%6d %9.5f %9.6f %6d %9.6f %9.6f %6d %s %s %-6s  ",  \
						Nidx, Qatm[i,j], radi, Catm[i,j], Satm[i,j], Eatm[i,j], \
						Ntot, Name[i]"-"n+1"."j+1, \
						resName[resID[i,j]], Tatm[i, j]
						if(Ntot in ndxPro) print "Pro"
						if(Ntot in ndxLig) print "Lig"
					}
				}
			}
		}
	}

	function getRadi(tag, radType) { # mBondi from AMBER20/parmed/tools/changeradii.py
		radBondi["C" ]= 1.7  ; radBondi["H" ]= 1.2
		radBondi["N" ]= 1.55 ; radBondi["HC"]= 1.3
		radBondi["O" ]= 1.5  ; radBondi["HN"]= 1.3
		radBondi["F" ]= 1.5  ; radBondi["HP"]= 1.3
		radBondi["SI"]= 2.1  ; radBondi["HO"]= 0.8
		radBondi["P" ]= 1.85 ; radBondi["HS"]= 0.8
		radBondi["S" ]= 1.8  # radBondi["BR"]= 1.85
		radBondi["CL"]= 1.7  # radBondi["I" ]= 1.98

		tag=toupper(tag)
		if(length(tag)>=2) {
			if(radBondi[substr(tag,1,2)]) return radBondi[substr(tag,1,2)]
			else return radBondi[substr(tag,1,1)]
		} else {
			if(radBondi[tag]) return radBondi[tag]
			else return 1.5
		}
	}
'

[[ ! -f "$qrv" ]] && { echo -e "!!! ERROR !!! gmx dump Failed ! Check $err to find out why.\n"; exit; }

if [[ $isCAS -eq 1 ]]; then
################################################################################
# CAS: qrv pdb
################################################################################

awk -v qrv=$qrv -v cas="$cas" -v RS="\r?\n" '
	BEGIN {
		txt="NME Z ALA A ARG R ASN N ASP D ASH D CYS C CYM C CYX C "
		txt=txt"GLU E GLH E GLN Q GLY G HID H HIE H HIP H ILE I LEU L LYS K LYN K "
		txt=txt"MET M PHE F PRO P SER S THR T TRP W TYR Y VAL V"
		n=split(txt, arr)
		for(i=1; i<=n/2; i++) resnm[arr[2*i-1]]=arr[2*i]

		split(cas,arr,",")
		for(i in arr) {
			if(index(arr[i],":")) {
				n=split(arr[i], arr2, ":")
				s=1; if(n>2) s=arr2[3]
				for(j=arr2[1]; j<=arr2[2]; j+=s) resCAS[j]=1
			} else resCAS[arr[i]]=1
		}

		getline < qrv
		getline Atyp < qrv
		for(i=0; i<Atyp; i++) getline < qrv
		while(getline < qrv) {
			res=$(NF-2); gsub(/[0-9]+/,"",res)
			if(res=="ALA") { n=$(NF-1);
				Qala[n]=$2; Rala[n]=$3; Cala[n]=$4
				Sala[n]=$5; Eala[n]=$6
			}
		}
		close(qrv)

		out=qrv; sub(".qrv", "_CAS.qrv", out)
		getline < qrv;      print     >out
		getline Atyp < qrv; print Atyp>out
		for(i=0; i<Atyp; i++) { getline < qrv; print >out }
		while(getline < qrv) {
			res=$(NF-2); gsub(/[0-9]+/,"",res)
			idx=$(NF-2); gsub(/^0*/,"",idx); gsub(/[A-Z]+/,"",idx)
			atm=$(NF-1)

			if(!resCAS[idx]) print >out
			else {
				if(res=="GLY") {
				# CH3
				} else if(res=="PRO") {
				# HN
				} else if(  ( res=="SER" &&  atm=="OG"                ) \
						 || ( res=="THR" && (atm=="OG1"||atm=="CG2") )  \
						 || ((res=="ILE"||res=="VAL") && (atm=="CG1"||atm=="CG2") ) \
						 || ((res=="CYS"||res=="CYM"||res=="CYX") && atm=="SG"    ) \
						 || ( atm=="CG") ) atm="HB3"

				if(atm=="N"  || atm=="H"   || atm=="HN"  \
				|| atm=="C"  || atm=="O"                 \
				|| atm=="CA" || atm=="CB"  || atm=="HA"  || atm=="HCA" \
				|| atm=="HB" || atm=="HB1" || atm=="HB2" || atm=="HB3") {
					$2=Qala[atm]; $5=Sala[atm]
					$3=Rala[atm]; $6=Eala[atm]
					$4=Cala[atm]; $9=sprintf("%05d%s", idx, resnm[res]"2"resnm["ALA"])
					printf "%6d %9.5f %9.6f %6d %9.6f %9.6f %6d %s %s %-6s  %s\n",
						$1, $2, $3, $4, $5, $6, $7, $8, $9, atm, $11 >out
				}
			}
		}
		close(qrv)
	}

	/^ATOM/ {
		ATOM=substr($0,1,6)
		INDX=substr($0,7,5)+0; isCAS[INDX]=0
		NAME=substr($0,13,4); name=NAME; gsub(/ /,"", name)
		RES =substr($0,18,3)
		CHN =substr($0,22,1);  if(CHN=" ") CHN="A"
		NUM =substr($0,23,4)+0
		X=substr($0,31,8); X += 0
		Y=substr($0,39,8); Y += 0
		Z=substr($0,47,8); Z += 0

		if(resCAS[NUM]) {
			Rhb=1.1
			if(name=="CB") { Xcb=X; Ycb=Y; Zcb=Z }

			if(RES=="GLY") {
			# CH3
			} else if(RES=="PRO") {
			# HN
			} else if(  ( RES=="SER" &&  name=="OG"                ) \
					 || ( RES=="THR" && (name=="OG1"||name=="CG2") ) \
					 || ((RES=="ILE"||RES=="VAL") && (name=="CG1"||name=="CG2") ) \
					 || ((RES=="CYS"||RES=="CYM"||RES=="CYX") && name=="SG"     ) \
					 || ( name=="CG") ) {
				NAME=" HB "; name="HB"
				dx=X-Xcb; dy=Y-Ycb; dz=Z-Zcb
				r=sqrt(dx*dx+dy*dy+dz*dz)
				X=Xcb+Rhb*dx/r
				Y=Ycb+Rhb*dy/r
				Z=Zcb+Rhb*dz/r
			}

			RES=resnm[RES]"2"resnm["ALA"]
			if(!(name=="N"  || name=="H"  || name=="HN"  \
			  || name=="CA" || name=="HA"                \
			  || name=="CB" || name=="HB" || name=="HB1" || name=="HB2" || name=="HB3" \
			  || name=="C"  || name=="O")) next
		}
		printf("%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f\n", \
			ATOM, INDX, NAME, RES, CHN, NUM, X, Y, Z)
		next
	}
	{print}
' $pdb > ${pdb/.pdb/_CAS.pdb}

qrv=${qrv/.qrv/_CAS.qrv}
pdb=${pdb/.pdb/_CAS.pdb}

fi

echo -e ">> 2. generate qrv file: OK !\n"

fi; if [[ $step -le 3 ]]; then
################################################################################
# 3. MM-PBSA计算: pdb->pqr, 输出apbs, 计算MM, APBS
# 3. run MM-PBSA, pdb->pqr, apbs, then calculate MM, PB, SA
################################################################################

echo -e ">> 3. run MM-PBSA calculatons"
dt=$(awk '/t=/{n++;sub(/.*t=/,"");sub(/step=.*/,"");t[n]=$0;if(n==2){print t[n]-t[1];exit}}' $pdb)
awk -v pid=_$pid          -v qrv=$qrv           -v dt="$dt"     \
	-v apbs="$apbs"       -v useDH=$useDH       -v useTS=$useTS \
	-v PBEset="$PBEset"   -v PBAset="$PBAset"                   \
	-v meshType=$meshType -v gridType=$gridType -v gmem=$gmem   \
	-v fadd=$fadd         -v cfac=$cfac         -v df=$df       \
	-v withLig=$withLig   -v RS="\r?\n" '
	BEGIN {
		getline < qrv
		getline Atyp < qrv
		for(i=0; i<Atyp; i++) {
			getline < qrv
			for(j=0; j<Atyp; j++) { C6[i, j]=$(2+2*j); C12[i,j]=$(3+2*j) }
		}
		n=0
		while(getline < qrv) { n++
			Qatm[$1]=$2; Ratm[$1]=$3; Catm[$1]=$4
			if($NF=="Pro") { Npro++; if(Npro==1) Ipro=n
				ndxPro[$1]++; resPro[Npro]="P~"$(NF-2)
			}
			if($NF=="Lig") { Nlig++; if(Nlig==1) Ilig=n
				ndxLig[$1]++; resLig[Nlig]="L~"$(NF-2)
			}
		}
		close(qrv)
		Ncom=Npro+Nlig

		txt=PBEset;
		gsub(/#[^\n]*\n/, "\n", txt);
		split(txt, arr, "\n")
		Nion=0
		for(i in arr) {
			split(arr[i], arr2)
			if(arr2[1]~/temp/) temp=arr2[2]
			if(arr2[1]~/pdie/) pdie=arr2[2]
			if(arr2[1]~/sdie/) sdie=arr2[2]
			if(arr2[1]~/ion/)  { Nion++; Qion[Nion]=arr2[3]; Cion[Nion]=arr2[5] }
		}

		Iion=0
		for(i=1; i<=Nion; i++) Iion += Cion[i]*Qion[i]^2
		eps0=8.854187812800001e-12
		kb=1.380649e-23
		Na=6.02214076e+23
		qe=1.602176634e-19
		RT2kJ=8.314462618*temp/1E3
		kap=1E-10/sqrt(eps0*kb*temp*sdie/(Iion*qe^2*Na*1E3))

		PBEset0=PBEset; sub(/sdie +[0-9]*\.*[0-9]*/, "sdie  1", PBEset0)

		txt=PBAset; split(txt, arr, "\n")
		for(i in arr) {
			if(arr[i]~/gamma +const/) {
				sub(/.*gamma +con[^ ]+ +/, "", arr[i]);
				split(arr[i], arr2)
				gamma=arr2[1]; const=arr2[2]
			}
		}

		MAXPOS=1E9
		minX= MAXPOS; maxX=-MAXPOS;
		minY= MAXPOS; maxY=-MAXPOS;
		minZ= MAXPOS; maxZ=-MAXPOS

		fmt=sprintf("%.9f",dt/1E3)
		sub(/0*$/,"",fmt);sub(/.*\./,"",fmt)
		fmt="~%."length(fmt)"fns"
	}

	/REMARK/ {next}
	/TITLE/ {
		if(Nfrm) {
			close(Fname[Nfrm]"_com.pqr")
			close(Fname[Nfrm]"_pro.pqr")
			if(withLig) close(Fname[Nfrm]"_lig.pqr")
		}
		Fout=FILENAME
		txt=$0; sub(/.*t= */,"",txt); sub(/ .*/,"",txt)
		txt=sprintf(fmt, txt/1E3);
		sub(".pdb", txt, Fout)
		Nfrm++
		Fname[Nfrm]=Fout

		minXpro[Nfrm]= MAXPOS; minXlig[Nfrm]= MAXPOS;
		minYpro[Nfrm]= MAXPOS; minYlig[Nfrm]= MAXPOS;
		minZpro[Nfrm]= MAXPOS; minZlig[Nfrm]= MAXPOS

		maxXpro[Nfrm]=-MAXPOS; maxXlig[Nfrm]=-MAXPOS
		maxYpro[Nfrm]=-MAXPOS; maxYlig[Nfrm]=-MAXPOS
		maxZpro[Nfrm]=-MAXPOS; maxZlig[Nfrm]=-MAXPOS
	}
	/^ATOM/ {
		ATOM=substr($0,1,6)
		INDX=substr($0,7,5)+0
		NAME=substr($0,13,4)
		RES =substr($0,18,3)
		CHN =substr($0,22,1); if(CHN=" ") CHN="A"
		NUM =substr($0,23,4)
		X   =substr($0,31,8); X += 0
		Y   =substr($0,39,8); Y += 0
		Z   =substr($0,47,8); Z += 0
		r=Ratm[INDX]

		txt=sprintf("%-6s%5d %-4s %3s %s%4d    %8.3f %8.3f %8.3f %12.6f %12.6f", \
			ATOM, INDX, NAME, RES, CHN, NUM, X, Y, Z, Qatm[INDX], r)

		if(INDX in ndxPro) {
			print txt > Fout"_pro.pqr"
			minXpro[Nfrm]=min(minXpro[Nfrm], X-r); maxXpro[Nfrm]=max(maxXpro[Nfrm], X+r)
			minYpro[Nfrm]=min(minYpro[Nfrm], Y-r); maxYpro[Nfrm]=max(maxYpro[Nfrm], Y+r)
			minZpro[Nfrm]=min(minZpro[Nfrm], Z-r); maxZpro[Nfrm]=max(maxZpro[Nfrm], Z+r)
		}

		if(withLig) {
			print txt > Fout"_com.pqr"
			if(INDX in ndxLig) {
				print txt > Fout"_lig.pqr"
				minXlig[Nfrm]=min(minXlig[Nfrm], X-r); maxXlig[Nfrm]=max(maxXlig[Nfrm], X+r)
				minYlig[Nfrm]=min(minYlig[Nfrm], Y-r); maxYlig[Nfrm]=max(maxYlig[Nfrm], Y+r)
				minZlig[Nfrm]=min(minZlig[Nfrm], Z-r); maxZlig[Nfrm]=max(maxZlig[Nfrm], Z+r)
			}
		}

		minXcom[Nfrm]=min(minXpro[Nfrm], minXlig[Nfrm]); maxXcom[Nfrm]=max(maxXpro[Nfrm], maxXlig[Nfrm])
		minYcom[Nfrm]=min(minYpro[Nfrm], minYlig[Nfrm]); maxYcom[Nfrm]=max(maxYpro[Nfrm], maxYlig[Nfrm])
		minZcom[Nfrm]=min(minZpro[Nfrm], minZlig[Nfrm]); maxZcom[Nfrm]=max(maxZpro[Nfrm], maxZlig[Nfrm])

		minX=min(minX, minXcom[Nfrm]); maxX=max(maxX, maxXcom[Nfrm])
		minY=min(minY, minYcom[Nfrm]); maxY=max(maxY, maxYcom[Nfrm])
		minZ=min(minZ, minZcom[Nfrm]); maxZ=max(maxZ, maxZcom[Nfrm])
		next
	}

	END {
		close(Fname[Nfrm]"_com.pqr")
		close(Fname[Nfrm]"_pro.pqr")
		if(withLig) close(Fname[Nfrm]"_lig.pqr")

		kJcou=1389.35457520287
		Rcut=1E10              # large enough

		for(i=1; i<=Npro; i++) dE[resPro[i]]=0
		for(i=1; i<=Nlig; i++) dE[resLig[i]]=0
		Nres=asorti(dE, Tres)

		txt="#Frame   "; txtDH=txt
		for(i=1; i<=Nres; i++) {
			ii=Tres[i]; sub(/~0+/, "~", ii)
			txt   = txt""sprintf("%10s", ii)
			txtDH = txtDH""sprintf("%21s", ii"(with DH)")
		}
		print txt > pid"~resPBSA.dat"
		print txt > pid"~resPBSA_PB.dat"
		print txt > pid"~resPBSA_SA.dat"
		if(withLig) {
			print txt   > pid"~resMM_VDW.dat"
			print txtDH > pid"~resMM.dat"
			print txtDH > pid"~resMM_COU.dat"
			print txtDH > pid"~res_MMPBSA.dat"
		}

		print "   #Frame      Binding( with DH ) "         \
			 "|    MM    ( with DH )    PB        SA     " \
			 "|   COU    ( with DH )     VDW   "            \
			 "|       PBcom        PBpro        PBlig  "   \
			 "|    SAcom     SApro     SAlig"    >> pid"~MMPBSA.dat"

		maxstr=0
		for(fr=1; fr<=Nfrm; fr++) maxstr=max(maxstr, length(Fname[fr]))

		if(withLig) {
			for(i=1; i<=Nres; i++) {
				tag=Tres[i]
				 dGres[tag,0]=0;  dGres[tag,1]=0
				 dHres[tag,0]=0;  dHres[tag,1]=0
				 MMres[tag,0]=0;  MMres[tag,1]=0
				COUres[tag,0]=0; COUres[tag,1]=0
				VDWres[tag]=0;   dPBres[tag]=0;  dSAres[tag]=0
			}
		}

		for(fr=1; fr<=Nfrm; fr++) {
			t_start = systime(); fflush("")

			Fout=Fname[fr]
			printf "   >> Frame %3d/%d: %-"maxstr"s", fr, Nfrm, Fout

			txt=Fout"_pro.pqr"; if(withLig) txt=Fout"_com.pqr";
			close(txt)
			n=0;
			while(getline < txt) { n++;
				qr[n]=$2; type[n]=$3; res[n]=$4;
				x[n]=$(NF-4);    y[n]=$(NF-3);   z[n]=$(NF-2)
				resID[n]=$(NF-5); gsub(/[A-Z]+/, "", resID[n])
			}
			close(txt)

			# MM
			if(withLig) {
				for(i=1; i<=Npro; i++) { dEcou[resPro[i]]=0; dEcouDH[resPro[i]]=0; dEvdw[resPro[i]]=0 }
				for(i=1; i<=Nlig; i++) { dEcou[resLig[i]]=0; dEcouDH[resLig[i]]=0; dEvdw[resLig[i]]=0 }
				for(i=1; i<=Npro; i++) {
					ii=i+Ipro-1
					qi=Qatm[qr[ii]]; ci=Catm[qr[ii]]
					xi=x[ii]; yi=y[ii]; zi=z[ii]
					for(j=1; j<=Nlig; j++) {
						jj=j+Ilig-1; cj=Catm[qr[jj]]
						r=sqrt( (xi-x[jj])^2+(yi-y[jj])^2+(zi-z[jj])^2 )
						if(r<Rcut) {
							t=1/(.1*r)^6
							Ecou   = qi*Qatm[qr[jj]]/r;
							EcouDH = Ecou*exp(-kap*r)
							Evdw = (C12[ci,cj]*t-C6[ci,cj])*t
							# print i,j,C12[ci,cj],C6[ci,cj],r,t,Evdw,xi,yi,zi,x[jj],y[jj],z[jj]> "Evdw.txt"
							dEcou[  resPro[i]] += Ecou;   dEcou[  resLig[j]] += Ecou
							dEcouDH[resPro[i]] += EcouDH; dEcouDH[resLig[j]] += EcouDH
							dEvdw[  resPro[i]] += Evdw;   dEvdw[  resLig[j]] += Evdw
							print i,j,qi,Qatm[qr[jj]],Ecou,dEcou[  resPro[i]]> "Ecou.txt"
						}
					}
				}

				Ecou=0; Evdw=0
				for(i in dEcou) {
					# print i,dEcou[  i]> "Ecou1.txt"
					dEcou[  i] *= kJcou/(2*pdie); Ecou   += dEcou[i];
					dEcouDH[i] *= kJcou/(2*pdie); EcouDH += dEcouDH[i];
					dEvdw[  i] /= 2;              Evdw   += dEvdw[i]
				}
			}

			# PBSA
			if(withLig) print "read\n" \
				"  mol pqr "Fout"_com.pqr\n" \
				"  mol pqr "Fout"_pro.pqr\n" \
				"  mol pqr "Fout"_lig.pqr\n" \
				"end\n\n" > Fout".apbs"
			else        print "read\n" \
				"  mol pqr "Fout"_pro.pqr\n" \
				"end\n\n" > Fout".apbs"

			if(meshType==0) { # GMXPBSA
				if(withLig) print \
					dimAPBS(Fout"_com", 1, minX, maxX, minY, maxY, minZ, maxZ), \
					dimAPBS(Fout"_pro", 2, minX, maxX, minY, maxY, minZ, maxZ), \
					dimAPBS(Fout"_lig", 3, minX, maxX, minY, maxY, minZ, maxZ)  > Fout".apbs"
				else        print \
					dimAPBS(Fout"_pro", 1, minX, maxX, minY, maxY, minZ, maxZ)  > Fout".apbs"
				print minX,minY,minZ,maxX,maxY,maxZ >> "XYZ.txt"
			} else if(meshType==1) { # g_mmpbsa
				if(withLig) print \
					dimAPBS(Fout"_com", 1, minXcom[fr], maxXcom[fr], minYcom[fr], maxYcom[fr], minZcom[fr], maxZcom[fr]), \
					dimAPBS(Fout"_pro", 2, minXpro[fr], maxXpro[fr], minYpro[fr], maxYpro[fr], minZpro[fr], maxZpro[fr]), \
					dimAPBS(Fout"_lig", 3, minXlig[fr], maxXlig[fr], minYlig[fr], maxYlig[fr], minZlig[fr], maxZlig[fr])  > Fout".apbs"
				else       print \
					dimAPBS(Fout"_pro", 1, minXpro[fr], maxXpro[fr], minYpro[fr], maxYpro[fr], minZpro[fr], maxZpro[fr])  > Fout".apbs"
			}
			close(Fout".apbs")

			cmd=apbs" "Fout".apbs > "Fout".out 2>&1";
			cmd | getline err; close(cmd)

			txt=Fout".out";
			while(getline < txt ) {
				if(index($0, "CALCULATION #")) {
					if(index($0, "("Fout"_com")) { t=1; n=Ncom }
					if(index($0, "("Fout"_pro")) { t=2; n=Npro }
					if(index($0, "("Fout"_lig")) { t=3; n=Nlig }
					if(index($0, "~VAC)")) t += 10
					if(index($0, "~SAS)")) t += 20
					while(getline < txt) {
						if(t<20 && index($0, "Per-atom energies:") \
						|| t>20 && index($0, "Solvent Accessible Surface Area")) break
					}

					for(i=1; i<=n; i++) {
						getline <txt;
						if(t<20) r=$3; else r=$NF
						if(t<10)       Esol[t%10, i]=r
						else if(t<20)  Evac[t%10, i]=r
						else if(t<30)  Esas[t%10, i]=gamma*r+const/n
					}
				}
			}
			close(txt)

			PBcom=0; SAcom=0;
			PBpro=0; SApro=0;
			PBlig=0; SAlig=0;
			for(i=1; i<=Ncom; i++) { Esol[1,i] -= Evac[1,i]; PBcom += Esol[1,i]; SAcom += Esas[1,i] }
			for(i=1; i<=Npro; i++) { Esol[2,i] -= Evac[2,i]; PBpro += Esol[2,i]; SApro += Esas[2,i] }
			for(i=1; i<=Nlig; i++) { Esol[3,i] -= Evac[3,i]; PBlig += Esol[3,i]; SAlig += Esas[3,i] }

			for(i=1; i<=Npro; i++) { dPBres[resPro[i]]=0; dSAres[resPro[i]]=0 }
			for(i=1; i<=Nlig; i++) { dPBres[resLig[i]]=0; dSAres[resLig[i]]=0 }
			for(i=1; i<=Npro; i++) {
				dPBres[resPro[i]] += Esol[1, Ipro+i-1]-Esol[2, i]
				dSAres[resPro[i]] += Esas[1, Ipro+i-1]-Esas[2, i]
			}
			for(i=1; i<=Nlig; i++) {
				dPBres[resLig[i]] += Esol[1, Ilig+i-1]-Esol[3, i]
				dSAres[resLig[i]] += Esas[1, Ilig+i-1]-Esas[3, i]
			}

			preK=-1; if(withLig) preK=1
			vdw[fr]=Evdw
			pb[fr]=preK*(PBcom-PBpro-PBlig)
			sa[fr]=preK*(SAcom-SApro-SAlig)
			cou[fr,0]=Ecou;
			cou[fr,1]=EcouDH
			mm[fr,0]=cou[fr,0]+vdw[fr];
			mm[fr,1]=cou[fr,1]+vdw[fr]
			dh[fr,0]=preK*mm[fr,0] + pb[fr]+sa[fr]
			dh[fr,1]=preK*mm[fr,1] + pb[fr]+sa[fr]
			printf "%-12s %9.3f(%9.3f) | %9.3f(%9.3f) %9.3f %9.3f " \
				  "| %9.3f(%9.3f) %9.3f | %12.3f %12.3f %12.3f | %9.3f %9.3f %9.3f\n", \
				Fout, dh[fr,0], dh[fr,1], mm[fr,0], mm[fr,1], pb[fr], sa[fr], \
				Ecou, EcouDH, Evdw, PBcom, PBpro, PBlig, SAcom, SApro, SAlig >> pid"~MMPBSA.dat"

			total_time = (Nfrm-fr)*(systime() - t_start)
			printf "  MM-PBSA(with DH) = %9.3f(%9.3f) kJ/mol    Estimated Time Remaining: ~ %0dh:%02dm:%02ds\n", \
				dh[fr,0], dh[fr,1], int(total_time/3600), int(total_time%3600/60), total_time%60

			fmt ="%s%9.3f %s"
			fmt2="%s%9.3f(%9.3f) %s"
			for(i=1; i<=Nres; i++) {
				ii="";  if(i==1) ii=sprintf("%-9s", Fout)
				txt=""; if(i==Nres) txt="\n"
				tag=Tres[i]
				if(withLig) {
					printf fmt,  ii, dEvdw[tag], txt               > pid"~resMM_VDW.dat"
					printf fmt2, ii, dEcou[tag], dEcouDH[tag], txt > pid"~resMM_COU.dat"
					printf fmt2, ii, dEcou[tag]+dEvdw[tag], dEcouDH[tag]+dEvdw[tag], txt  > pid"~resMM.dat"
					printf fmt2, ii, dEcou[tag]+dEvdw[tag]+dPBres[tag]+dSAres[tag],
									 dEcouDH[tag]+dEvdw[tag]+dPBres[tag]+dSAres[tag], txt > pid"~res_MMPBSA.dat"

					COUres[tag,0] += dEcou[tag]
					COUres[tag,1] += dEcouDH[tag]
					VDWres[tag]+= dEvdw[tag];
				}
				PBres[tag] += dPBres[tag]
				SAres[tag] += dSAres[tag]
				printf fmt, ii, preK*(dPBres[tag]+dSAres[tag]), txt > pid"~resPBSA.dat"
				printf fmt, ii, preK*(dPBres[tag]), txt             > pid"~resPBSA_PB.dat"
				printf fmt, ii, preK*(dSAres[tag]), txt             > pid"~resPBSA_SA.dat"
			}

			fmt="%s%6s%6s\n"
			for(i=1; i<=Npro; i++) {
				ii=Ipro+i-1
				txt=sprintf("%-6s%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f", \
					"ATOM", ii, type[ii], res[ii], resID[ii], x[ii], y[ii], z[ii])

				tag=resPro[i]
				printf fmt, txt, fixfmt(preK*dPBres[tag],6), fixfmt(preK*dPBres[tag],6) > Fout"~resPBSA_PB.pdb"
				printf fmt, txt, fixfmt(preK*dSAres[tag],6), fixfmt(preK*dSAres[tag],6) > Fout"~resPBSA_SA.pdb"
				if(withLig) {
					printf fmt, txt, fixfmt(dEcou[tag],6),              fixfmt(dEcouDH[tag],6)            > Fout"~resMM_COU.pdb"
					printf fmt, txt, fixfmt(dEvdw[tag],6),              fixfmt(dEvdw[tag],6)              > Fout"~resMM_VDW.pdb"
					printf fmt, txt, fixfmt(dEcou[tag]+dEvdw[tag],6),   fixfmt(dEcouDH[tag]+dEvdw[tag],6) > Fout"~resMM.pdb"
					printf fmt, txt, fixfmt(dPBres[tag]+dSAres[tag],6), fixfmt(dPBres[tag]+dSAres[tag],6) > Fout"~resPBSA.pdb"
					printf fmt, txt, fixfmt(dEcou[tag]+dEvdw[tag]+preK*(dPBres[tag]+dSAres[tag]),6),
									 fixfmt(dEcouDH[tag]+dEvdw[tag]+preK*(dPBres[tag]+dSAres[tag]),6)     > Fout"~res_MMPBSA.pdb"
				} else
					printf fmt, txt, fixfmt(preK*(dPBres[tag]+dSAres[tag]),6),
									 fixfmt(preK*(dPBres[tag]+dSAres[tag]),6) > Fout"~resPBSA.pdb"
			}
			for(i=1; i<=Nlig; i++) {
				ii=Ilig+i-1
				txt=sprintf("%-6s%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f", \
					 "ATOM", ii, type[ii], res[ii], resID[ii], x[ii], y[ii], z[ii])
				tag=resLig[i]
				printf fmt, txt, fixfmt(dEcou[tag]+dEvdw[tag]+preK*(dPBres[tag]+dSAres[tag]),6),
								 fixfmt(dEcouDH[tag]+dEvdw[tag]+preK*(dPBres[tag]+dSAres[tag]),6)     > Fout"~res_MMPBSA.pdb"
				printf fmt, txt, fixfmt(preK*dPBres[tag],6),        fixfmt(preK*dPBres[tag],6)        > Fout"~resPBSA_PB.pdb"
				printf fmt, txt, fixfmt(preK*dSAres[tag],6),        fixfmt(preK*dSAres[tag],6)        > Fout"~resPBSA_SA.pdb"
				printf fmt, txt, fixfmt(dEcou[tag],6),              fixfmt(dEcouDH[tag],6)            > Fout"~resMM_COU.pdb"
				printf fmt, txt, fixfmt(dEvdw[tag],6),              fixfmt(dEvdw[tag],6)              > Fout"~resMM_VDW.pdb"
				printf fmt, txt, fixfmt(dEcou[tag]+dEvdw[tag],6),   fixfmt(dEcouDH[tag]+dEvdw[tag],6) > Fout"~resMM.pdb"
				printf fmt, txt, fixfmt(dPBres[tag]+dSAres[tag],6), fixfmt(dPBres[tag]+dSAres[tag],6) > Fout"~resPBSA.pdb"
			}

			close(Fout"~resMM.pdb")
			close(Fout"~resMM_COU.pdb")
			close(Fout"~resMM_VDW.pdb")
			close(Fout"~resPBSA.pdb")
			close(Fout"~resPBSA_PB.pdb")
			close(Fout"~resPBSA_SA.pdb")
			close(Fout"~res_MMPBSA.pdb")
		}

		print "#mol res       MM-PBSA(with DH  )        MM(with DH  )      PBSA       " \
			  "COU(with DH  )       VDW        PB        SA" > pid"~res.dat"
		fmt="%-12s %9.3f(%9.3f) %9.3f(%9.3f) %9.3f %9.3f(%9.3f) %9.3f %9.3f %9.3f\n"
		for(i=1; i<=Nres; i++) {
			tag=Tres[i]

			VDWres[tag]   /= Nfrm
			PBres[tag]    /= Nfrm
			SAres[tag]    /= Nfrm
			COUres[tag,0] /= Nfrm
			COUres[tag,1] /= Nfrm

			MMres[tag,0]=COUres[tag,0]+VDWres[tag]
			MMres[tag,1]=COUres[tag,1]+VDWres[tag]
			dHres[tag,0]=MMres[tag,0]+PBres[tag]+SAres[tag]
			dHres[tag,1]=MMres[tag,1]+PBres[tag]+SAres[tag]

			txt=tag; sub(/L~0+/, "Lig ", txt); sub(/P~0+/, "Pro ", txt)
			printf fmt, txt,
				dHres[tag,0], dHres[tag,1],
				MMres[tag,0], MMres[tag,1],
				PBres[tag]+SAres[tag],
				COUres[tag,0], COUres[tag,1], VDWres[tag], \
				PBres[tag], SAres[tag]  > pid"~res.dat"

			ii="";  if(i==1) ii=sprintf("---------------------------------------\n%-9s", "mean")
			txt=""; if(i==Nres) txt="\n"
			if(withLig) {
				printf "%s%9.3f(%9.3f) %s", ii, dHres[tag,0], dHres[tag,1]  , txt > pid"~res_MMPBSA.dat"
				printf "%s%9.3f(%9.3f) %s", ii, MMres[tag,0], MMres[tag,1]  , txt > pid"~resMM.dat"
				printf "%s%9.3f %s",        ii, VDWres[tag]                 , txt > pid"~resMM_VDW.dat"
				printf "%s%9.3f(%9.3f) %s", ii, COUres[tag,0], COUres[tag,1], txt > pid"~resMM_COU.dat"
			}
			printf "%s%9.3f %s",        ii, preK*(PBres[tag]+SAres[tag]), txt > pid"~resPBSA.dat"
			printf "%s%9.3f %s",        ii, preK*(PBres[tag])           , txt > pid"~resPBSA_PB.dat"
			printf "%s%9.3f %s",        ii, preK*(SAres[tag])           , txt > pid"~resPBSA_SA.dat"
		}

		fmt="%s%6s%6s\n"
		for(i=1; i<=Npro; i++) {
			ii=Ipro+i-1
			txt=sprintf("%-6s%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f", \
				"ATOM", ii, type[ii], res[ii], resID[ii], x[ii], y[ii], z[ii])

			tag=resPro[i]
			if(withLig) {
				printf fmt, txt, fixfmt(dHres[tag,0],6),             fixfmt(dHres[tag,1],6)                 > pid"~res_MMPBSA.pdb"
				printf fmt, txt, fixfmt(MMres[tag,0]+VDWres[tag],6), fixfmt(MMres[tag,1]+VDWres[tag],6)     > pid"~resMM.pdb"
				printf fmt, txt, fixfmt(COUres[tag,0],6),            fixfmt(COUres[tag,1],6)                > pid"~resMM_COU.pdb"
				printf fmt, txt, fixfmt(VDWres[tag],6),              fixfmt(VDWres[tag],6)                  > pid"~resMM_VDW.pdb"
			}
			printf fmt, txt, fixfmt(preK*(PBres[tag]+SAres[tag]),6), fixfmt(preK*(PBres[tag]+SAres[tag]),6) > pid"~resPBSA.pdb"
			printf fmt, txt, fixfmt(preK*PBres[tag]+SAres[tag],6),   fixfmt(preK*PBres[tag]+SAres[tag],6)   > pid"~resPBSA_PB.pdb"
			printf fmt, txt, fixfmt(preK*SAres[tag],6),              fixfmt(preK*SAres[tag],6)              > pid"~resPBSA_SA.pdb"
		}
		for(i=1; i<=Nlig; i++) {
			ii=Ilig+i-1
			txt=sprintf("%-6s%5d %-4s %3s A%4d    %8.3f%8.3f%8.3f", \
				 "ATOM", ii, type[ii], res[ii], resID[ii], x[ii], y[ii], z[ii])
			tag=resLig[i]
			printf fmt, txt, fixfmt(dHres[tag,0],6),                 fixfmt(dHres[tag,1],6)                 > pid"~res_MMPBSA.pdb"
			printf fmt, txt, fixfmt(MMres[tag,0]+VDWres[tag],6),     fixfmt(MMres[tag,1]+VDWres[tag],6)     > pid"~resMM.pdb"
			printf fmt, txt, fixfmt(COUres[tag,0],6),                fixfmt(COUres[tag,1],6)                > pid"~resMM_COU.pdb"
			printf fmt, txt, fixfmt(VDWres[tag],6),                  fixfmt(VDWres[tag],6)                  > pid"~resMM_VDW.pdb"
			printf fmt, txt, fixfmt(preK*(PBres[tag]+SAres[tag]),6), fixfmt(preK*(PBres[tag]+SAres[tag]),6) > pid"~resPBSA.pdb"
			printf fmt, txt, fixfmt(preK*PBres[tag]+SAres[tag],6),   fixfmt(preK*PBres[tag]+SAres[tag],6)   > pid"~resPBSA_PB.pdb"
			printf fmt, txt, fixfmt(preK*SAres[tag],6),              fixfmt(preK*SAres[tag],6)              > pid"~resPBSA_SA.pdb"
		}

		 dH[0]=0;  dH[1]=0
		 MM[0]=0;  MM[1]=0
		COU[0]=0; COU[1]=0
		VDW=0;    PB=0;     SA=0
		for(i=1; i<=Nfrm; i++) {
			dH[0] += dh[i,0]/Nfrm;  dH[1] += dh[i,1]/Nfrm
			MM[0] += mm[i,0]/Nfrm;  MM[1] += mm[i,1]/Nfrm
			COU[0]+= cou[i,0]/Nfrm; COU[1]+= cou[i,1]/Nfrm
			VDW += vdw[i]/Nfrm
			PB  += pb[i]/Nfrm
			SA  += sa[i]/Nfrm
		}

		TdS[0]=0; TdS[1]=0 # 1 for DH
		if(withLig) {
			for(i=1; i<=Nfrm; i++) {
				TdS[0] += exp((mm[i,0]-MM[0])/RT2kJ)/Nfrm
				TdS[1] += exp((mm[i,1]-MM[1])/RT2kJ)/Nfrm
			}
		}
		TdS[0]=-RT2kJ*log(TdS[0]); TdS[1]=-RT2kJ*log(TdS[1]);
		 dG[0]=dH[0]-TdS[0];        dG[1]=dH[1]-TdS[1];
		 Ki[0]=exp(dG[0]/RT2kJ);    Ki[1]=exp(dG[1]/RT2kJ)

		printf "----------------------------------" \
			"|------------------------------------------" \
			"|--------------------------------|\n" \
			"%-12s %9.3f(%9.3f) | %9.3f(%9.3f) %9.3f %9.3f | %9.3f(%9.3f) %9.3f |\n"   \
			"%-12s %9.3f(%9.3f) |\n----------------------------------|\n" \
			"%-12s %9.3f(%9.3f) kJ/mol = %9.3f(%9.3f) kcal/mol\n" \
			"%-12s %9.3E(%9.3E)  uM    = %9.3E(%9.3E)    nM\n" ,
			"  dH=",   dH[0],     dH[1],     MM[0], MM[1], PB, SA, COU[0], COU[1], VDW,
			"-TdS=", -TdS[0],   -TdS[1],
			"  dG=",   dG[0],     dG[1],     dG[0]/4.184, dG[1]/4.184,
			"  Ki=",   Ki[0]*1E6, Ki[1]*1E6, Ki[0]*1E9,   Ki[1]*1E9     >> pid"~MMPBSA.dat"

		print "+==============================================================+\n"\
			"             dG = RTln(Ki) = -RTln(KA)                          \n" \
			"             Ki = 1/KA = exp(dG/RT) = IC50 = EC50               \n" \
			"+==============================================================+\n" \
			"| M(mol/L) |  m/u/pM | dG(kcal/mol) | dG(kJ/mol) | Affinity    |\n" \
			"|==============================================================|\n" \
			"| 0.1      | 100  mM |    -1.364    |   -5.708   |             |\n" \
			"| 0.01     |  10  mM |    -2.728    |  -11.416   |             |\n" \
			"| 0.001    |   1  mM |    -4.093    |  -17.124   |             |\n" \
			"| 0.0001   | 100  uM |    -5.457    |  -22.832   | Weak        |\n" \
			"|--------------------------------------------------------------|\n" \
			"| 0.00001  |  10  uM |    -6.821    |  -28.540   |             |\n" \
			"| 1.0E-06  |   1  uM |    -8.185    |  -34.248   | Medium      |\n" \
			"|--------------------------------------------------------------|\n" \
			"| 1.0E-07  | 100  nM |    -9.550    |  -39.956   |             |\n" \
			"| 1.0E-08  |  10  nM |   -10.914    |  -45.664   |             |\n" \
			"| 1.0E-09  |   1  nM |   -12.278    |  -51.372   | Strong      |\n" \
			"|--------------------------------------------------------------|\n" \
			"| 1.0E-10  | 100  pM |   -13.642    |  -57.080   |             |\n" \
			"| 1.0E-11  |  10  pM |   -15.007    |  -62.788   |             |\n" \
			"| 1.0E-12  |   1  pM |   -16.371    |  -68.496   | Very strong |\n" \
			"+==============================================================+\n" >> pid"~MMPBSA.dat"
	}

	function dimAPBS(file, Imol, minX, maxX, minY, maxY, minZ, maxZ) {

		lenX=max(maxX-minX, 0.1); cntX=(maxX+minX)/2
		lenY=max(maxY-minY, 0.1); cntY=(maxY+minY)/2
		lenZ=max(maxZ-minZ, 0.1); cntZ=(maxZ+minZ)/2
		cX  =lenX*cfac;           fX  =min(cX, lenX+fadd)
		cY  =lenY*cfac;           fY  =min(cY, lenY+fadd)
		cZ  =lenZ*cfac;           fZ  =min(cZ, lenZ+fadd)

		levN=4    # 划分级别
		t=2^(levN+1)
		nX=round(fX/df)-1; nX=max(t*round(nX/t)+1, 33)
		nY=round(fY/df)-1; nY=max(t*round(nY/t)+1, 33)
		nZ=round(fZ/df)-1; nZ=max(t*round(nZ/t)+1, 33)

		if(gridType==0) { # GMXPBSA method
			fpre=1; cfac=1.7
			fX=lenX+2*fadd; cX=fX*cfac; nX=t*(int(fX/(t*df))+1+fpre)+1
			fY=lenY+2*fadd; cY=fY*cfac; nY=t*(int(fY/(t*df))+1+fpre)+1
			fZ=lenZ+2*fadd; cZ=fZ*cfac; nZ=t*(int(fZ/(t*df))+1+fpre)+1
		}

		MGset="mg-auto"
		mem = 200*nX*nY*nZ/1024./1024. # MB

#		npX=nX; npY=nY; npZ=nZ
#		gmem=4000
#		ofrac=0.1
#		if(mem>=gmem) {
#			while(mem>gmem) {
#				maxN=max(npX, max(npY, npZ))
#					 if(maxN==npX) npX = t*((npX-1)/t-1)+1
#				else if(maxN==npY) npY = t*((npY-1)/t-1)+1
#				else if(maxN==npZ) npZ = t*((npZ-1)/t-1)+1
#				mem = 200*npX*npY*npZ/1024./1024
#			}

#			t=nX/npX; if(t>1) npX = int(t*(1+2*ofrac) + 1.0);
#			t=nY/npY; if(t>1) npY = int(t*(1+2*ofrac) + 1.0);
#			t=nZ/npZ; if(t>1) npZ = int(t*(1+2*ofrac) + 1.0);
#			MGset="mg-para\n  ofrac "ofrac"\n  pdime "npX" "npY" "npZ
#		}

		XYZset="  "MGset \
			"\n  mol "Imol \
			"\n  dime   "nX"  "nY"  "nZ"        # 格点数目, 所需内存: "mem" MB"  \
			"\n  cglen  "cX"  "cY"  "cZ"        # 粗略格点长度" \
			"\n  fglen  "fX"  "fY"  "fZ"        # 细密格点长度" \
			"\n  fgcent "cntX"  "cntY"  "cntZ"  # 细密格点中心" \
			"\n  cgcent "cntX"  "cntY"  "cntZ"  # 粗略格点中心"

		return \
			"ELEC name "file"\n" \
			XYZset "\n" \
			PBEset "\n" \
			"end\n\n" \
			"ELEC name "file"~VAC\n" \
			XYZset  "\n" \
			PBEset0 "\n" \
			"end\n\n" \
			"APOLAR name "file"~SAS\n" \
			"  mol "Imol"\n" \
			PBAset"\n" \
			"end\n\n" \
			"print elecEnergy "file" - "file"~VAC end\n" \
			"print apolEnergy "file"~SAS end\n\n"
	}
	function min(x, y)    { return x<y ? x : y }
	function max(x, y)    { return x>y ? x : y }
	function round(x)     { return int(x+0.5)  }
	function fixfmt(x,m, nint,npre,nexp,str){
		str=sprintf("%.9f", x)
		sub(/^0+/,   "", str);
		sub(/^-0+/,  "-", str);
		nint=index(str, ".")-1
		if(-10^m<x && x<-10^(-m) || 10^(-m)<x &&x<10^m ) {
			str=sprintf("%."(m-nint-1)"f", x)
			sub(/^0+/,   "", str);
			sub(/^-0+/,  "-", str);
		} else {
			str=sprintf("%.9E", x)
			npre=index(str, ".")
			sub(/.*E\+0*$/, "",  str);
			sub(/.*E\+0*/, "E",  str);
			sub(/.*E-0*/,  "E-", str)
			nexp=length(str)
			str=sprintf("%"(m-nexp)"."(m-nexp-npre)"E", x)
			sub(/E\+0*$/, "",  str);
			sub(/E\+0*/, "E",  str);
			sub(/E-0*/,  "E-", str)
		}
		return str
	}
' $pdb
echo -e "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
fi

################################################################################
# 4. 删除临时文件
# 4. remove intermediate files
################################################################################
# rm -f io.mc \#*
#rm -f _$pid*