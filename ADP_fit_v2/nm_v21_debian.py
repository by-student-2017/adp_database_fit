from scipy import optimize
import numpy
import numpy as np
import commands
import sys
#----------------------------------------------------------------------
fix_atom_position = 0 # 0: Off, 1: On
#----------------------------------------------------------------------
file_tmp = 'ADP_code_v21.tmp'
file_inp = 'ADP_code_v21'

cif2cell_adress = "cif2cell"

commands.getoutput("setenv OMP_NUM_THREADS 1")
num_core = commands.getoutput("grep 'core id' /proc/cpuinfo | sort -u | wc -l")
#lammps_adress = "mpirun -np "+str(num_core)+" --allow-run-as-root lmp"
#pwscf_adress = "mpirun -np "+str(num_core)+" --allow-run-as-root pw.x"
lammps_adress = "mpirun -np "+str(num_core)+" lmp"
pwscf_adress = "mpirun -np "+str(num_core)+" pw.x"
#lammps_adress = "mpirun -np 2 lmp"
#pwscf_adress = "mpirun -np 2 pw.x"

satom = commands.getoutput("grep \"atomtype\" ADP.input | sed -e \"s/.*=//\" -e \"s/'//g\"")

commands.getoutput("chmod +x ./cfg2vasp/cfg2vasp")
commands.getoutput("chmod +x pwscf2force")
commands.getoutput("cp in.lmp_temp_v2 in.lmp_temp")
commands.getoutput("chmod +x setinp")
commands.getoutput("./setinp")
commands.getoutput("mkdir cfg")
commands.getoutput("mkdir work")
commands.getoutput("echo -n > energy.dat")

natom = 5000
fxl = numpy.ones(int(natom)+1)
fyl = numpy.ones(int(natom)+1)
fzl = numpy.ones(int(natom)+1)
fxp = numpy.ones(int(natom)+1)
fyp = numpy.ones(int(natom)+1)
fzp = numpy.ones(int(natom)+1)

print "use struct.dat"
struct = commands.getoutput("awk '{if($1==\""+str(satom)+"\"){print $0}}' struct.dat")
struct_list = struct.split()
ntemp = int((len(struct_list)-1)/3 - 1)
temp = []
stru = []
weig = []
#if float(struct_list[3*ntemp+1]) <= 1073.0 :
#  ntemp = ntemp + 1
#  struct_list.append(1273.0)
#  struct_list.append("L")
#  struct_list.append(1.0)
for i in range(ntemp+1):
  temp.append(float(struct_list[3*i+1]))
  stru.append(struct_list[3*i+2])
  weig.append(float(struct_list[3*i+3]))
  t = temp[i]
  s = stru[i]
  commands.getoutput("cp in.lmp in.lmp_"+str(t)+"K")
  commands.getoutput("sed -i 's/YYYY/"+str(t)+"/' in.lmp_"+str(t)+"K")
  commands.getoutput("cp ./data/data.in."+str(s)+" data.in_"+str(t)+"K")
print "temperature: ",temp
print "structure  : ",stru
print "weight     : ",weig
#----------------------------------------------------------------------
print "read parameters from ADP_code_v21.init"
nline = commands.getoutput("grep -n "+str(satom)+" ADP_code_v21.init | head -1 | sed -e \"s/:.*//g\"")
print "read line: "+nline
check_satom = commands.getoutput("awk '{if(NR=="+str(nline)+"+0){print $1}}' ADP_code_v21.init | head -1")
print "fit element: "+check_satom
# fitting parameters
x0  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+1){print $1}}' ADP_code_v21.init | head -1"))
x1  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+2){print $1}}' ADP_code_v21.init | head -1"))
x2  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+3){print $1}}' ADP_code_v21.init | head -1"))
x3  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+4){print $1}}' ADP_code_v21.init | head -1"))
x4  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+5){print $1}}' ADP_code_v21.init | head -1"))
x5  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+6){print $1}}' ADP_code_v21.init | head -1"))
x6  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+7){print $1}}' ADP_code_v21.init | head -1"))
x7  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+8){print $1}}' ADP_code_v21.init | head -1"))
x8  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+9){print $1}}' ADP_code_v21.init | head -1"))
x9  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+10){print $1}}' ADP_code_v21.init | head -1"))
x10 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+11){print $1}}' ADP_code_v21.init | head -1"))
x11 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+12){print $1}}' ADP_code_v21.init | head -1"))
x12 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+13){print $1}}' ADP_code_v21.init | head -1"))
x13 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+14){print $1}}' ADP_code_v21.init | head -1"))
x14 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+15){print $1}}' ADP_code_v21.init | head -1"))
x15 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+16){print $1}}' ADP_code_v21.init | head -1"))
x16 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+17){print $1}}' ADP_code_v21.init | head -1"))
x17 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+18){print $1}}' ADP_code_v21.init | head -1"))
x18 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+19){print $1}}' ADP_code_v21.init | head -1"))
x19 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+20){print $1}}' ADP_code_v21.init | head -1"))
x20 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+21){print $1}}' ADP_code_v21.init | head -1"))
x21 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+22){print $1}}' ADP_code_v21.init | head -1"))
x22 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+23){print $1}}' ADP_code_v21.init | head -1"))
x23 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+24){print $1}}' ADP_code_v21.init | head -1"))
x24 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+25){print $1}}' ADP_code_v21.init | head -1"))
x25 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+26){print $1}}' ADP_code_v21.init | head -1"))
x26 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+27){print $1}}' ADP_code_v21.init | head -1"))
x27 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+30){print $1}}' ADP_code_v21.init | head -1"))
x28 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+31){print $1}}' ADP_code_v21.init | head -1"))
x29 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+32){print $1}}' ADP_code_v21.init | head -1"))
# u
x30 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+34){print $1}}' ADP_code_v21.init | head -1"))
x31 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+35){print $1}}' ADP_code_v21.init | head -1"))
x32 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+36){print $1}}' ADP_code_v21.init | head -1"))
x33 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+37){print $1}}' ADP_code_v21.init | head -1"))
x34 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+38){print $1}}' ADP_code_v21.init | head -1"))
x35 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+39){print $1}}' ADP_code_v21.init | head -1"))
x36 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+40){print $1}}' ADP_code_v21.init | head -1"))
x37 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+41){print $1}}' ADP_code_v21.init | head -1"))
x38 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+42){print $1}}' ADP_code_v21.init | head -1"))
x39 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+43){print $1}}' ADP_code_v21.init | head -1"))
x40 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+44){print $1}}' ADP_code_v21.init | head -1"))
x41 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+45){print $1}}' ADP_code_v21.init | head -1"))
x42 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+46){print $1}}' ADP_code_v21.init | head -1"))
x43 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+47){print $1}}' ADP_code_v21.init | head -1"))
x44 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+48){print $1}}' ADP_code_v21.init | head -1"))
x45 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+49){print $1}}' ADP_code_v21.init | head -1"))
x46 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+50){print $1}}' ADP_code_v21.init | head -1"))
x47 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+51){print $1}}' ADP_code_v21.init | head -1"))
x48 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+52){print $1}}' ADP_code_v21.init | head -1"))
# w
x49 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+53){print $1}}' ADP_code_v21.init | head -1"))
x50 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+54){print $1}}' ADP_code_v21.init | head -1"))
x51 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+55){print $1}}' ADP_code_v21.init | head -1"))
x52 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+56){print $1}}' ADP_code_v21.init | head -1"))
x53 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+57){print $1}}' ADP_code_v21.init | head -1"))
x54 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+58){print $1}}' ADP_code_v21.init | head -1"))
x55 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+59){print $1}}' ADP_code_v21.init | head -1"))
x56 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+60){print $1}}' ADP_code_v21.init | head -1"))
x57 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+61){print $1}}' ADP_code_v21.init | head -1"))
x58 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+62){print $1}}' ADP_code_v21.init | head -1"))
x59 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+63){print $1}}' ADP_code_v21.init | head -1"))
x60 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+64){print $1}}' ADP_code_v21.init | head -1"))
x61 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+65){print $1}}' ADP_code_v21.init | head -1"))
x62 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+66){print $1}}' ADP_code_v21.init | head -1"))
x63 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+67){print $1}}' ADP_code_v21.init | head -1"))
x64 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+68){print $1}}' ADP_code_v21.init | head -1"))
x65 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+69){print $1}}' ADP_code_v21.init | head -1"))
x66 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+70){print $1}}' ADP_code_v21.init | head -1"))
x67 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+71){print $1}}' ADP_code_v21.init | head -1"))

z0 = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,
      x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48]
#      x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64,x65,x66,x67]
print "initial parameters: ",z0

count = 0
#----------------------------------------------------------------------
def f(z):
  
  print "------------------------"
  global count
  count += 1
  print count

  fi = open(file_tmp,'r')
  text = fi.read().replace('re',str(z[0]).replace("[","").replace("]",""))
  text = text.replace('p0',str(z[1]).replace("[","").replace("]",""))
  text = text.replace('Frhoe',str(z[2]).replace("[","").replace("]",""))
  text = text.replace('pk1',str(z[3]).replace("[","").replace("]",""))
  text = text.replace('pk2',str(z[4]).replace("[","").replace("]",""))
  text = text.replace('pk3',str(z[5]).replace("[","").replace("]",""))
  text = text.replace('pk4',str(z[6]).replace("[","").replace("]",""))
  text = text.replace('pk5',str(z[7]).replace("[","").replace("]",""))
  text = text.replace('pk6',str(z[8]).replace("[","").replace("]",""))
  text = text.replace('pk7',str(z[9]).replace("[","").replace("]",""))
  text = text.replace('pk8',str(z[10]).replace("[","").replace("]",""))
  text = text.replace('Fi0',str(z[11]).replace("[","").replace("]",""))
  text = text.replace('Fi1',str(z[12]).replace("[","").replace("]",""))
  text = text.replace('Fi2',str(z[13]).replace("[","").replace("]",""))
  text = text.replace('Fi3',str(z[14]).replace("[","").replace("]",""))
  text = text.replace('Fi4',str(z[15]).replace("[","").replace("]",""))
  text = text.replace('Fi5',str(z[16]).replace("[","").replace("]",""))
  text = text.replace('Fi6',str(z[17]).replace("[","").replace("]",""))
  text = text.replace('Fm2',str(z[18]).replace("[","").replace("]",""))
  text = text.replace('Fm3',str(z[19]).replace("[","").replace("]",""))
  text = text.replace('Fm4',str(z[20]).replace("[","").replace("]",""))
  text = text.replace('Fm5',str(z[21]).replace("[","").replace("]",""))
  text = text.replace('Fm6',str(z[22]).replace("[","").replace("]",""))
  text = text.replace('Fn0',str(z[23]).replace("[","").replace("]",""))
  text = text.replace('Fn1',str(z[24]).replace("[","").replace("]",""))
  text = text.replace('Fn2',str(z[25]).replace("[","").replace("]",""))
  text = text.replace('Fn3',str(z[26]).replace("[","").replace("]",""))
  text = text.replace('p1',str(z[27]).replace("[","").replace("]",""))
  text = text.replace('he',str(z[28]).replace("[","").replace("]",""))
  text = text.replace('Frhol',str(z[29]).replace("[","").replace("]",""))
  # u
  text = text.replace('urhoe',str(z[30]).replace("[","").replace("]",""))
  text = text.replace('ui0',str(z[31]).replace("[","").replace("]",""))
  text = text.replace('ui1',str(z[32]).replace("[","").replace("]",""))
  text = text.replace('ui2',str(z[33]).replace("[","").replace("]",""))
  text = text.replace('ui3',str(z[34]).replace("[","").replace("]",""))
  text = text.replace('ui4',str(z[35]).replace("[","").replace("]",""))
  text = text.replace('ui5',str(z[36]).replace("[","").replace("]",""))
  text = text.replace('ui6',str(z[37]).replace("[","").replace("]",""))
  text = text.replace('um2',str(z[38]).replace("[","").replace("]",""))
  text = text.replace('um3',str(z[39]).replace("[","").replace("]",""))
  text = text.replace('um4',str(z[40]).replace("[","").replace("]",""))
  text = text.replace('um5',str(z[41]).replace("[","").replace("]",""))
  text = text.replace('um6',str(z[42]).replace("[","").replace("]",""))
  text = text.replace('un0',str(z[43]).replace("[","").replace("]",""))
  text = text.replace('un1',str(z[44]).replace("[","").replace("]",""))
  text = text.replace('un2',str(z[45]).replace("[","").replace("]",""))
  text = text.replace('un3',str(z[46]).replace("[","").replace("]",""))
  text = text.replace('urhol',str(z[47]).replace("[","").replace("]",""))
  text = text.replace('urhoh',str(z[48]).replace("[","").replace("]",""))
  # w
  text = text.replace('wrhoe',str(x49).replace("[","").replace("]",""))
  text = text.replace('wi0',str(x50).replace("[","").replace("]",""))
  text = text.replace('wi1',str(x51).replace("[","").replace("]",""))
  text = text.replace('wi2',str(x52).replace("[","").replace("]",""))
  text = text.replace('wi3',str(x53).replace("[","").replace("]",""))
  text = text.replace('wi4',str(x54).replace("[","").replace("]",""))
  text = text.replace('wi5',str(x55).replace("[","").replace("]",""))
  text = text.replace('wi6',str(x56).replace("[","").replace("]",""))
  text = text.replace('wm2',str(x57).replace("[","").replace("]",""))
  text = text.replace('wm3',str(x58).replace("[","").replace("]",""))
  text = text.replace('wm4',str(x59).replace("[","").replace("]",""))
  text = text.replace('wm5',str(x60).replace("[","").replace("]",""))
  text = text.replace('wm6',str(x61).replace("[","").replace("]",""))
  text = text.replace('wn0',str(x62).replace("[","").replace("]",""))
  text = text.replace('wn1',str(x63).replace("[","").replace("]",""))
  text = text.replace('wn2',str(x64).replace("[","").replace("]",""))
  text = text.replace('wn3',str(x65).replace("[","").replace("]",""))
  text = text.replace('wrhol',str(x66).replace("[","").replace("]",""))
  text = text.replace('wrhoh',str(x67).replace("[","").replace("]",""))
  fi.close

  with open(file_inp,'w') as f:
    print >> f, text

  commands.getoutput("./Zhou04_ADP_v21 < ADP.input")

  tdiffea = 0.0
  tdiffp  = 0.0
  tdifff  = 0.0
  for t in temp:
    print "---------------"
    print "Temperature: "+str(t)+" [K]"
    if count > 20000 or count % int(10000*2.718**(-count/10000)+1) == 1:
      commands.getoutput("mv data.in_"+str(t)+"K data.in")
      natom = commands.getoutput("awk '{if($2==\"atoms\"){print $1}}' data.in")
      commands.getoutput(lammps_adress+" < in.lmp_"+str(t)+"K")
      if (fix_atom_position == 1):
        commands.getoutput("cp ./cfg/run.0.cfg run.50.cfg")
      else:
        commands.getoutput("cp ./cfg/run.50.cfg run.50.cfg")
      commands.getoutput("./cfg2vasp/cfg2vasp run.50.cfg")
      commands.getoutput("python ./vasp2cif/vasp2cif.py run.50.vasp")
      commands.getoutput(cif2cell_adress+" run.50.vasp.cif --no-reduce -p pwscf --pwscf-pseudo-PSLibrary-libdr=\"./potentials\" --setup-all --k-resolution=0.48 --pwscf-force=yes --pwscf-stress=yes --pwscf-run-type=scf -o pw.in") 
      commands.getoutput("sed -i 's/\'pw\'/\'pw_"+str(t)+"K\'/g' pw.scf.in")
      commands.getoutput(pwscf_adress+" < pw.scf.in > pw.out")
      commands.getoutput(cif2cell_adress+" run.50.vasp.cif --no-reduce -p pwscf --pwscf-pseudo-PSLibrary-libdr=\"./potentials\" --setup-all --k-resolution=0.20 --pwscf-force=yes --pwscf-stress=yes --pwscf-run-type=scf -o pw.in") 
      commands.getoutput("sed -i 's/\'pw\'/\'pw_"+str(t)+"K\'/g' pw.scf.in")
      commands.getoutput(pwscf_adress+" < pw.scf.in > pw.out")
      commands.getoutput("./pwscf2force >> config_potfit_"+str(satom))
      commands.getoutput(cif2cell_adress+" run.50.vasp.cif --no-reduce -p lammps -o data_fix.in_"+str(t)+"K")
      commands.getoutput("cp data_fix.in_"+str(t)+"K data_fix.in")
      commands.getoutput(lammps_adress+" < in.lmp_fix")
      commands.getoutput("mv data.in.restart data.in_"+str(t)+"K")
    #
      commands.getoutput("./pwscf2force > config_"+str(t)+"K")
    else:
      commands.getoutput("cp data_fix.in_"+str(t)+"K data_fix.in")
      natom = commands.getoutput("awk '{if($2==\"atoms\"){print $1}}' data_fix.in")
      commands.getoutput(lammps_adress+" < in.lmp_fix")
    print "number of atoms: "+str(natom)

    # 1 bar = 0.0001 GPa
    # stress = -pressure
    #pxxl = commands.getoutput("awk '{if($1==\"pxxl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    #pyyl = commands.getoutput("awk '{if($1==\"pyyl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    #pzzl = commands.getoutput("awk '{if($1==\"pzzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    #pxyl = commands.getoutput("awk '{if($1==\"pxyl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    #pxzl = commands.getoutput("awk '{if($1==\"pxzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    #pyzl = commands.getoutput("awk '{if($1==\"pyzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    #pxxp = commands.getoutput("awk '{if($1==\"#S\"){print -$2}}' config")
    #pyyp = commands.getoutput("awk '{if($1==\"#S\"){print -$3}}' config")
    #pzzp = commands.getoutput("awk '{if($1==\"#S\"){print -$4}}' config")
    #pxyp = commands.getoutput("awk '{if($1==\"#S\"){print -$5}}' config")
    #pxzp = commands.getoutput("awk '{if($1==\"#S\"){print -$6}}' config")
    #pyzp = commands.getoutput("awk '{if($1==\"#S\"){print -$7}}' config")
    #diffpxx = (float(pxxl) - float(pxxp))/(float(pxxp)+0.000000101)*100.0/6.0
    #diffpyy = (float(pyyl) - float(pyyp))/(float(pyyp)+0.000000101)*100.0/6.0
    #diffpzz = (float(pzzl) - float(pzzp))/(float(pzzp)+0.000000101)*100.0/6.0
    #diffpxy = (float(pxyl) - float(pxyp))/(float(pxyp)+0.000000101)*100.0/6.0
    #diffpxz = (float(pxzl) - float(pxzp))/(float(pxzp)+0.000000101)*100.0/6.0
    #diffpyz = (float(pyzl) - float(pyzp))/(float(pyzp)+0.000000101)*100.0/6.0
    #diffp = abs(diffpxx) + abs(diffpyy) + abs(diffpzz) + abs(diffpxy) + abs(diffpxz) + abs(diffpyz)
    #print "lammps: "+str(pxxl)+", "+str(pyyl)+", "+str(pzzl)+", "+str(pxyl)+", "+str(pxzl)+", "+str(pyzl)+" [eV/A^3]"
    #print "pwscf:  "+str(pxxp)+", "+str(pyyp)+", "+str(pzzp)+", "+str(pxyp)+", "+str(pxzp)+", "+str(pyzp)+" [eV/A^3]"
    #print "P diff (%): "+str(diffp)
    #print "---------------"
    diffp = 0.0

    # force
    difffx = 0.0
    difffy = 0.0
    difffz = 0.0
    difff  = 0.0
    for i in range(int(natom)):
      fxl[i] = commands.getoutput("awk '{if(NR==10+"+str(i)+"){printf \"%10.8f\",$7}}' trajectory.lammpstrj")
      fyl[i] = commands.getoutput("awk '{if(NR==10+"+str(i)+"){printf \"%10.8f\",$8}}' trajectory.lammpstrj")
      fzl[i] = commands.getoutput("awk '{if(NR==10+"+str(i)+"){printf \"%10.8f\",$9}}' trajectory.lammpstrj")
      fxp[i] = commands.getoutput("awk '{if(NR==11+"+str(i)+"){print $5}}' config_"+str(t)+"K")
      fyp[i] = commands.getoutput("awk '{if(NR==11+"+str(i)+"){print $6}}' config_"+str(t)+"K")
      fzp[i] = commands.getoutput("awk '{if(NR==11+"+str(i)+"){print $7}}' config_"+str(t)+"K")
      difffx = (float(fxl[i]) - float(fxp[i]))/(float(fxp[i])+0.000000101)*100.0/3.0/float(natom)
      difffy = (float(fyl[i]) - float(fyp[i]))/(float(fyp[i])+0.000000101)*100.0/3.0/float(natom)
      difffz = (float(fzl[i]) - float(fzp[i]))/(float(fzp[i])+0.000000101)*100.0/3.0/float(natom)
      difff  = difff + abs(difffx) + abs(difffy) + abs(difffz)
    print "lammps: "+str(fxl[0])+" : "+str(fyl[0])+" : "+str(fzl[0])+" [eV/A]"
    print "PWscf: "+str(fxp[0])+" : "+str(fyp[0])+" : "+str(fzp[0])+" [eV/A]"
    print "force diff (%): "+str(difff)
    print "---------------"

    lammps_get_data = "grep \"Total Energy\" log.lammps | tail -1 | awk '{printf \"%-20.10f\",$4}'"
    lmpe = commands.getoutput(lammps_get_data)

    pwe = commands.getoutput("awk '{if($1==\"#E\"){print $2}}' config_"+str(t)+"K")
    pwe = float(pwe) * float(natom)

    print "lammps: "+str(lmpe)+" [eV]"

    print "PWscf:  "+str(pwe)+" [eV]"

    diffe = float(pwe) - float(lmpe)
    print "diff: "+str(diffe)+" [eV]"
    diffea = float(diffe)/float(natom)
    print "diff/atom: "+str(diffea)+" [eV/atom]"
    commands.getoutput("echo "+str(count)+" "+str(diffe)+" >> energy.dat")

    for itw in range(ntemp+1):
      if t == temp[itw]:
        wt = weig[itw]

    tdiffea = tdiffea + float(diffea)*float(wt)
    tdiffp  = tdiffp  + float(diffp)*float(wt)
    tdifff  = tdifff  + float(difff)*float(wt)

  diffb  = commands.getoutput("cat diff.dat")
  print "F, u and w boundary, diff: "+str(diffb)
  print "---------------"

  y = float(tdiffea)**2 + 1000*float(diffb)**2 + 0.0000002*abs(tdiffp)**2 + 0.0000010*abs(tdifff)**2


  print "Evaluate: ", y
  #print "Parameters: ", x
  print "Parameters: z0 = "+"[ "+str(z[0])+","+str(z[1])+","+str(z[2])+","+str(z[3])+","+str(z[4])+","+str(z[5])+","+str(z[6])+","+str(z[7])+","+str(z[8])+","+str(z[9])+","+str(z[10])+","+str(z[11])+","+str(z[12])+","+str(z[13])+","+str(z[14])+","+str(z[15])+","+str(z[16])+","+str(z[17])+","+str(z[18])+","+str(z[19])+","+str(z[20])+","+str(z[21])+","+str(z[22])+","+str(z[23])+","+str(z[24])+","+str(z[25])+","+str(z[26])+","+str(z[27])+","+str(z[28])+","+str(z[29])+","+str(z[30])+","+str(z[31])+","+str(z[32])+","+str(z[33])+","+str(z[34])+","+str(z[35])+","+str(z[36])+","+str(z[37])+","+str(z[38])+","+str(z[39])+","+str(z[40])+","+str(z[41])+","+str(z[42])+","+str(z[43])+","+str(z[44])+","+str(z[45])+","+str(z[46])+","+str(z[47])+","+str(z[48])+","+str(x[49])+","+str(x[50])+","+str(x[51])+","+str(x[52])+","+str(x[53])+","+str(x[54])+","+str(x[55])+","+str(x[56])+","+str(x[57])+","+str(x[58])+","+str(x[59])+","+str(x[60])+","+str(x[61])+","+str(x[62])+","+str(x[63])+","+str(x[64])+","+str(x[65])+","+str(x[66])+","+str(x[67])+" ]"
  print "------------------------"

  return y
#----------------------------------------------------------------------
res = optimize.minimize(f,z0,method='Nelder-Mead',options={'adaptive':True})
#res = optimize.minimize(f,x0,method='Nelder-Mead')
#res = optimize.minimize(f,x0,method='TNC')
#res = optimize.minimize(f,x0,method='Powell')
#res = optimize.minimize(f,x0,method='BFGS')
