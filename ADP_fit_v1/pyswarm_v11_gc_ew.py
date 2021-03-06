import random
#from scipy import optimize
#from fstpso import FuzzyPSO
from pyswarm import pso
import numpy
import numpy as np
import commands
import sys

print "Attention"
print "this code need Zhou04_ADP_EW"
print "  In the first time, need to command 'gfortran Zhou04_create_adp_v11.f -o Zhou04_ADP_v11'"
print "this code does not recommend element < Al"
print "all u parameters = 0.0"
print " "
#----------------------------------------------------------------------
fix_atom_position = 0 # 0: Off, 1: On
#----------------------------------------------------------------------
file_tmp = 'ADP_code_v11.tmp'
file_tmp_ew = 'ADP_code_ew_v11.tmp'
file_inp = 'ADP_code_v11'

cif2cell_adress = "cif2cell"

commands.getoutput("setenv OMP_NUM_THREADS 1")
num_core = commands.getoutput("grep 'core id' /proc/cpuinfo | sort -u | wc -l")
lammps_adress = "mpirun -np "+str(num_core)+" --allow-run-as-root lmp"
pwscf_adress = "mpirun -np "+str(num_core)+" --allow-run-as-root pw.x"
#lammps_adress = "mpirun -np "+str(num_core)+" lmp"
#pwscf_adress = "mpirun -np "+str(num_core)+" pw.x"
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
print "read parameters from ADP_code_v11.init"
nline = commands.getoutput("grep -n "+str(satom)+" ADP_code_v11.init | head -1 | sed -e \"s/:.*//g\"")
print "read line: "+nline
check_satom = commands.getoutput("awk '{if(NR=="+str(nline)+"+0){print $1}}' ADP_code_v11.init | head -1")
print "fit element: "+check_satom
# fitting parameters
x0  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+1){print $1}}' ADP_code_v11.init | head -1"))
x1  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+2){print $1}}' ADP_code_v11.init | head -1"))
x2  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+3){print $1}}' ADP_code_v11.init | head -1"))
x3  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+4){print $1}}' ADP_code_v11.init | head -1"))
x4  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+5){print $1}}' ADP_code_v11.init | head -1"))
x5  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+6){print $1}}' ADP_code_v11.init | head -1"))
x6  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+7){print $1}}' ADP_code_v11.init | head -1"))
x7  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+8){print $1}}' ADP_code_v11.init | head -1"))
x8  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+9){print $1}}' ADP_code_v11.init | head -1"))
x9  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+10){print $1}}' ADP_code_v11.init | head -1"))
x10 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+11){print $1}}' ADP_code_v11.init | head -1"))
x11 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+12){print $1}}' ADP_code_v11.init | head -1"))
x12 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+13){print $1}}' ADP_code_v11.init | head -1"))
x13 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+14){print $1}}' ADP_code_v11.init | head -1"))
x14 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+15){print $1}}' ADP_code_v11.init | head -1"))
x15 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+16){print $1}}' ADP_code_v11.init | head -1"))
x16 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+17){print $1}}' ADP_code_v11.init | head -1"))
x17 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+18){print $1}}' ADP_code_v11.init | head -1"))
x18 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+19){print $1}}' ADP_code_v11.init | head -1"))
x19 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+20){print $1}}' ADP_code_v11.init | head -1"))
x20 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+21){print $1}}' ADP_code_v11.init | head -1"))
x21 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+22){print $1}}' ADP_code_v11.init | head -1"))
x22 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+23){print $1}}' ADP_code_v11.init | head -1"))
x23 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+24){print $1}}' ADP_code_v11.init | head -1"))
x24 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+25){print $1}}' ADP_code_v11.init | head -1"))
x25 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+26){print $1}}' ADP_code_v11.init | head -1"))
x26 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+31){print $1}}' ADP_code_v11.init | head -1"))
# u
x27 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+33){print $1}}' ADP_code_v11.init | head -1"))
x28 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+34){print $1}}' ADP_code_v11.init | head -1"))
x29 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+35){print $1}}' ADP_code_v11.init | head -1"))
x30 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+36){print $1}}' ADP_code_v11.init | head -1"))
x31 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+37){print $1}}' ADP_code_v11.init | head -1"))
x32 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+38){print $1}}' ADP_code_v11.init | head -1"))
x33 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+39){print $1}}' ADP_code_v11.init | head -1"))
x34 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+40){print $1}}' ADP_code_v11.init | head -1"))
x35 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+41){print $1}}' ADP_code_v11.init | head -1"))
x36 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+42){print $1}}' ADP_code_v11.init | head -1"))
x37 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+43){print $1}}' ADP_code_v11.init | head -1"))
x38 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+44){print $1}}' ADP_code_v11.init | head -1"))
x39 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+45){print $1}}' ADP_code_v11.init | head -1"))
x40 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+46){print $1}}' ADP_code_v11.init | head -1"))
x41 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+47){print $1}}' ADP_code_v11.init | head -1"))
x42 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+48){print $1}}' ADP_code_v11.init | head -1"))
x43 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+49){print $1}}' ADP_code_v11.init | head -1"))
x44 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+50){print $1}}' ADP_code_v11.init | head -1"))
x45 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+51){print $1}}' ADP_code_v11.init | head -1"))
# w
x46 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+52){print $1}}' ADP_code_v11.init | head -1"))
x47 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+53){print $1}}' ADP_code_v11.init | head -1"))
x48 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+54){print $1}}' ADP_code_v11.init | head -1"))
x49 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+55){print $1}}' ADP_code_v11.init | head -1"))
x50 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+56){print $1}}' ADP_code_v11.init | head -1"))
x51 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+57){print $1}}' ADP_code_v11.init | head -1"))
x52 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+58){print $1}}' ADP_code_v11.init | head -1"))
x53 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+59){print $1}}' ADP_code_v11.init | head -1"))
x54 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+60){print $1}}' ADP_code_v11.init | head -1"))
x55 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+61){print $1}}' ADP_code_v11.init | head -1"))
x56 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+62){print $1}}' ADP_code_v11.init | head -1"))
x57 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+63){print $1}}' ADP_code_v11.init | head -1"))
x58 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+64){print $1}}' ADP_code_v11.init | head -1"))
x59 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+65){print $1}}' ADP_code_v11.init | head -1"))
x60 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+66){print $1}}' ADP_code_v11.init | head -1"))
x61 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+67){print $1}}' ADP_code_v11.init | head -1"))
x62 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+68){print $1}}' ADP_code_v11.init | head -1"))
x63 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+69){print $1}}' ADP_code_v11.init | head -1"))
x64 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+70){print $1}}' ADP_code_v11.init | head -1"))

# u
fi0 = open(file_tmp,'r')
text0 = fi0.read().replace('urhoe',str(x27).replace("[","").replace("]",""))
text0 = text0.replace('ui0',str(x28).replace("[","").replace("]",""))
text0 = text0.replace('ui1',str(x29).replace("[","").replace("]",""))
text0 = text0.replace('ui2',str(x30).replace("[","").replace("]",""))
text0 = text0.replace('ui3',str(x31).replace("[","").replace("]",""))
text0 = text0.replace('ui4',str(x32).replace("[","").replace("]",""))
text0 = text0.replace('ui5',str(x33).replace("[","").replace("]",""))
text0 = text0.replace('ui6',str(x34).replace("[","").replace("]",""))
text0 = text0.replace('um2',str(x35).replace("[","").replace("]",""))
text0 = text0.replace('um3',str(x36).replace("[","").replace("]",""))
text0 = text0.replace('um4',str(x37).replace("[","").replace("]",""))
text0 = text0.replace('um5',str(x38).replace("[","").replace("]",""))
text0 = text0.replace('um6',str(x39).replace("[","").replace("]",""))
text0 = text0.replace('un0',str(x40).replace("[","").replace("]",""))
text0 = text0.replace('un1',str(x41).replace("[","").replace("]",""))
text0 = text0.replace('un2',str(x42).replace("[","").replace("]",""))
text0 = text0.replace('un3',str(x43).replace("[","").replace("]",""))
text0 = text0.replace('urhol',str(x44).replace("[","").replace("]",""))
text0 = text0.replace('urhoh',str(x45).replace("[","").replace("]",""))
with open(file_tmp_ew,'w') as f:
  print >> f, text0

x = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,
      x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64]
print "initial parameters: ", x

count = 0
#----------------------------------------------------------------------
n_gene = int(len(x)) # number of parameters
min_ind = numpy.ones(n_gene) * -1.0
max_ind = numpy.ones(n_gene) *  1.0
min_max_ind = []
for i in range(n_gene):
  #min_ind[i] = b1[i][0]
  #max_ind[i] = b1[i][1]
  if (i < 27):
    min_ind[i] = float(x[i]) - float(abs(x[i]))*0.15
    max_ind[i] = float(x[i]) + float(abs(x[i]))*0.15 + 0.000001
  else:
    min_ind[i] = float(x[i]) - float(abs(x[i]))*0.07
    max_ind[i] = float(x[i]) + float(abs(x[i]))*0.07 + 0.000001
  #print "search area of paramter "+str(i)+": "+str(min_ind[i])+" | "+str(max_ind[i])
  min_max_ind.append([min_ind[i],max_ind[i]])
print "search area of paramter ", min_max_ind
#----------------------------------------------------------------------
def example_fitness( individual ):
  
  print "------------------------"
  global count
  count += 1
  print count

  fi = open(file_tmp_ew,'r')
  text = fi.read().replace('re',str(individual[0]).replace("[","").replace("]",""))
  text = text.replace('fe',str(individual[1]).replace("[","").replace("]",""))
  text = text.replace('Frhoe',str(individual[2]).replace("[","").replace("]",""))
  text = text.replace('Frhos',str(individual[3]).replace("[","").replace("]",""))
  text = text.replace('alpha',str(individual[4]).replace("[","").replace("]",""))
  text = text.replace('beta',str(individual[5]).replace("[","").replace("]",""))
  text = text.replace('Ap',str(individual[6]).replace("[","").replace("]",""))
  text = text.replace('Bp',str(individual[7]).replace("[","").replace("]",""))
  text = text.replace('cai',str(individual[8]).replace("[","").replace("]",""))
  text = text.replace('lambda',str(individual[9]).replace("[","").replace("]",""))
  text = text.replace('Fi0',str(individual[10]).replace("[","").replace("]",""))
  text = text.replace('Fi1',str(individual[11]).replace("[","").replace("]",""))
  text = text.replace('Fi2',str(individual[12]).replace("[","").replace("]",""))
  text = text.replace('Fi3',str(individual[13]).replace("[","").replace("]",""))
  text = text.replace('Fi4',str(individual[14]).replace("[","").replace("]",""))
  text = text.replace('Fi5',str(individual[15]).replace("[","").replace("]",""))
  text = text.replace('Fi6',str(individual[16]).replace("[","").replace("]",""))
  text = text.replace('Fm2',str(individual[17]).replace("[","").replace("]",""))
  text = text.replace('Fm3',str(individual[18]).replace("[","").replace("]",""))
  text = text.replace('Fm4',str(individual[19]).replace("[","").replace("]",""))
  text = text.replace('Fm5',str(individual[20]).replace("[","").replace("]",""))
  text = text.replace('Fm6',str(individual[21]).replace("[","").replace("]",""))
  text = text.replace('Fn0',str(individual[22]).replace("[","").replace("]",""))
  text = text.replace('Fn1',str(individual[23]).replace("[","").replace("]",""))
  text = text.replace('Fn2',str(individual[24]).replace("[","").replace("]",""))
  text = text.replace('Fn3',str(individual[25]).replace("[","").replace("]",""))
  text = text.replace('Frhol',str(individual[26]).replace("[","").replace("]",""))
  # w
  text = text.replace('wrhoe',str(individual[27]).replace("[","").replace("]",""))
  text = text.replace('wi0',str(individual[28]).replace("[","").replace("]",""))
  text = text.replace('wi1',str(individual[29]).replace("[","").replace("]",""))
  text = text.replace('wi2',str(individual[30]).replace("[","").replace("]",""))
  text = text.replace('wi3',str(individual[31]).replace("[","").replace("]",""))
  text = text.replace('wi4',str(individual[32]).replace("[","").replace("]",""))
  text = text.replace('wi5',str(individual[33]).replace("[","").replace("]",""))
  text = text.replace('wi6',str(individual[34]).replace("[","").replace("]",""))
  text = text.replace('wm2',str(individual[35]).replace("[","").replace("]",""))
  text = text.replace('wm3',str(individual[36]).replace("[","").replace("]",""))
  text = text.replace('wm4',str(individual[37]).replace("[","").replace("]",""))
  text = text.replace('wm5',str(individual[38]).replace("[","").replace("]",""))
  text = text.replace('wm6',str(individual[39]).replace("[","").replace("]",""))
  text = text.replace('wn0',str(individual[40]).replace("[","").replace("]",""))
  text = text.replace('wn1',str(individual[41]).replace("[","").replace("]",""))
  text = text.replace('wn2',str(individual[42]).replace("[","").replace("]",""))
  text = text.replace('wn3',str(individual[43]).replace("[","").replace("]",""))
  text = text.replace('wrhol',str(individual[44]).replace("[","").replace("]",""))
  text = text.replace('wrhoh',str(individual[45]).replace("[","").replace("]",""))
  fi.close

  with open(file_inp,'w') as f:
    print >> f, text

  commands.getoutput("./Zhou04_ADP_v11 < ADP.input")
  diffb  = commands.getoutput("cat diff.dat")
  if diffb == "nan" or abs(float(diffb)) >= 0.15/(1+float(count)/900):
    y = 99999.99999 + random.random()
    if count == 1:
      count -= 1
    print "skip this potential, because of bad boundary."
    return y

  tdiffea = 0.0
  tdiffp  = 0.0
  tdifff  = 0.0
  for t in temp:
    print "---------------"
    print "Temperature: "+str(t)+" [K]"
    if (count % 9000) == 1:
      commands.getoutput("mv data.in_"+str(t)+"K data.in")
      natom = commands.getoutput("awk '{if($2==\"atoms\"){print $1}}' data.in")
      commands.getoutput(lammps_adress+" < in.lmp_"+str(t)+"K")
      error_flag1 = ""
      error_flag2 = ""
      error_flag3 = ""
      error_flag1 = commands.getoutput("grep 'Total wall time' log.lammps")
      error_flag2 = commands.getoutput("grep 'nan' log.lammps")
      error_flag3 = commands.getoutput("grep 'ERROR' log.lammps")
      print error_flag1, error_flag2, error_flag3
      if error_flag1 != "" and error_flag2 == "" and error_flag3 == "":
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
        commands.getoutput("./data2cfg/lmp_data2cfg data.in "+satom)
        commands.getoutput("mv data.in.cfg run.0.cfg")
        commands.getoutput("./cfg2vasp/cfg2vasp run.0.cfg")
        commands.getoutput("python ./vasp2cif/vasp2cif.py run.0.vasp")
        commands.getoutput(cif2cell_adress+" run.0.vasp.cif --no-reduce -p pwscf --pwscf-pseudo-PSLibrary-libdr=\"./potentials\" --setup-all --k-resolution=0.48 --pwscf-force=yes --pwscf-stress=yes --pwscf-run-type=scf -o pw.in")
        commands.getoutput("sed -i 's/\'pw\'/\'pw_"+str(t)+"K\'/g' pw.scf.in")
        commands.getoutput(pwscf_adress+" < pw.scf.in > pw.out")
        commands.getoutput(cif2cell_adress+" run.0.vasp.cif --no-reduce -p pwscf --pwscf-pseudo-PSLibrary-libdr=\"./potentials\" --setup-all --k-resolution=0.20 --pwscf-force=yes --pwscf-stress=yes --pwscf-run-type=scf -o pw.in")
        commands.getoutput("sed -i 's/\'pw\'/\'pw_"+str(t)+"K\'/g' pw.scf.in")
        commands.getoutput(pwscf_adress+" < pw.scf.in > pw.out")
        commands.getoutput("./pwscf2force >> config_potfit_"+str(satom))
        commands.getoutput(cif2cell_adress+" run.0.vasp.cif --no-reduce -p lammps -o data_fix.in_"+str(t)+"K")
        commands.getoutput("cp data_fix.in_"+str(t)+"K data_fix.in")
        commands.getoutput("./pwscf2force > config_"+str(t)+"K")
    else:
      commands.getoutput("cp data_fix.in_"+str(t)+"K data_fix.in")
      natom = commands.getoutput("awk '{if($2==\"atoms\"){print $1}}' data_fix.in")
      commands.getoutput(lammps_adress+" < in.lmp_fix")
      error_flag1 = ""
      error_flag2 = ""
      error_flag3 = ""
      error_flag1 = commands.getoutput("grep 'Total wall time' log.lammps")
      error_flag2 = commands.getoutput("grep 'nan' log.lammps")
      error_flag3 = commands.getoutput("grep 'ERROR' log.lammps")
      print error_flag1, error_flag2, error_flag3
   
    if error_flag1 != "" and error_flag2 == "" and error_flag3 == "":

      print "number of atoms: "+str(natom)

      # 1 bar = 0.0001 GPa
      # stress = -pressure
      #pxxl = commands.getoutput("awk '{if($1==\"pxxl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
      #pyyl = commands.getoutput("awk '{if($1==\"pyyl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
      #pzzl = commands.getoutput("awk '{if($1==\"pzzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
      #pxyl = commands.getoutput("awk '{if($1==\"pxyl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
      #pxzl = commands.getoutput("awk '{if($1==\"pxzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
      #pyzl = commands.getoutput("awk '{if($1==\"pyzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
      #adddata = 0.0
      #if abs(float(pxxl)) <= 0.000000001 and abs(float(pyyl)) <= 0.000000001 and abs(float(pzzl)) <= 0.000000001 and abs(float(pxyl)) <= 0.000000001 and abs(float(pxzl)) <= 0.000000001 and abs(float(pyzl)) <= 0.000000001:
      #  pxxl = "99999.99999"
      #  pyyl = "99999.99999"
      #  pzzl = "99999.99999"
      #  pxyl = "99999.99999"
      #  pxzl = "99999.99999"
      #  pyzl = "99999.99999"
      #  adddata = random.random()
      #pxxp = commands.getoutput("awk '{if($1==\"#S\"){print $2}}' config_"+str(t)+"K")
      #pyyp = commands.getoutput("awk '{if($1==\"#S\"){print $3}}' config_"+str(t)+"K")
      #pzzp = commands.getoutput("awk '{if($1==\"#S\"){print $4}}' config_"+str(t)+"K")
      #pxyp = commands.getoutput("awk '{if($1==\"#S\"){print $5}}' config_"+str(t)+"K")
      #pxzp = commands.getoutput("awk '{if($1==\"#S\"){print $6}}' config_"+str(t)+"K")
      #pyzp = commands.getoutput("awk '{if($1==\"#S\"){print $7}}' config_"+str(t)+"K")
      #diffpxx = (float(pxxl) - float(pxxp))/(float(pxxp)+0.000000101)*100.0/6.0
      #diffpyy = (float(pyyl) - float(pyyp))/(float(pyyp)+0.000000101)*100.0/6.0
      #diffpzz = (float(pzzl) - float(pzzp))/(float(pzzp)+0.000000101)*100.0/6.0
      #diffpxy = (float(pxyl) - float(pxyp))/(float(pxyp)+0.000000101)*100.0/6.0
      #diffpxz = (float(pxzl) - float(pxzp))/(float(pxzp)+0.000000101)*100.0/6.0
      #diffpyz = (float(pyzl) - float(pyzp))/(float(pyzp)+0.000000101)*100.0/6.0
      #diffp = abs(diffpxx) + abs(diffpyy) + abs(diffpzz) + abs(diffpxy) + abs(diffpxz) + abs(diffpyz) + adddata
      #print "lammps: "+str(pxxl)+", "+str(pyyl)+", "+str(pzzl)+", "+str(pxyl)+", "+str(pxzl)+", "+str(pyzl)+" [eV/A^3]"
      #print "PWscf:  "+str(pxxp)+", "+str(pyyp)+", "+str(pzzp)+", "+str(pxyp)+", "+str(pxzp)+", "+str(pyzp)+" [eV/A^3]"
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
        if fxl[i] == 0.0 and fyl[i] == 0.0 and fzl[i] == 0.0:
          adddata = random.random()
          fxl[i] = 99999.99999 + adddata
          fyl[i] = 99999.99999 + adddata
          fzl[i] = 99999.99999 + adddata
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

    else:
      adddata = random.random()
      tdiffea = 99999.99999 + adddata
      tdiffp  = 99999.99999 + adddata
      tdifff  = 99999.99999 + adddata

  if error_flag1 != "" and error_flag2 == "" and error_flag3 == "":
    diffb  = commands.getoutput("cat diff.dat")
  else:
    diffb = 99999.99999 + random.random()
  print "F boundary, diff: "+str(diffb)
  print "---------------"

  y = float(tdiffea)**2 + 1000*float(diffb)**2 + 0.0000002*abs(tdiffp)**2 + 0.0000010*abs(tdifff)**2

  print "Evaluate: ", y
  #print "Parameters: ", x
  print "Parameters: x = "+"[ "+str(individual[0])+","+str(individual[1])+","+str(individual[2])+","+str(individual[3])+","+str(individual[4])+","+str(individual[5])+","+str(individual[6])+","+str(individual[7])+","+str(individual[8])+","+str(individual[9])+","+str(individual[10])+","+str(individual[11])+","+str(individual[12])+","+str(individual[13])+","+str(individual[14])+","+str(individual[15])+","+str(individual[16])+","+str(individual[17])+","+str(individual[18])+","+str(individual[19])+","+str(individual[20])+","+str(individual[21])+","+str(individual[22])+","+str(individual[23])+","+str(individual[24])+","+str(individual[25])+","+str(individual[26])+","+str(individual[27])+","+str(individual[28])+","+str(individual[29])+","+str(individual[30])+","+str(individual[31])+","+str(individual[32])+","+str(individual[33])+","+str(individual[34])+","+str(individual[35])+","+str(individual[36])+","+str(individual[37])+","+str(individual[38])+","+str(individual[39])+","+str(individual[40])+","+str(individual[41])+","+str(individual[42])+","+str(individual[43])+","+str(individual[44])+","+str(individual[45])+" ]"
  print "------------------------"

  return y
#----------------------------------------------------------------------
##PySwarm
print ('*'*65)
print ('Example minimization of 4th-order banana function (no constraints)')
#xopt1, fopt1 = pso(example_fitness, min_ind, max_ind)
xopt1, fopt1 = pso(example_fitness, min_ind, max_ind, maxiter=200)
#xopt1, fopt1 = pso(example_fitness, min_ind, max_ind, ieqcons=[],
#                   f_ieqcons=None, args=(), kwargs={},
#                   swarmsize=100, omega=0.5, phip=0.5, phig=0.5, maxiter=100,
#                   minstep=1e-8, minfunc=1e-8, debug=False)
# https://pythonhosted.org/pyswarm/
print ('The optimum is at:')
print ('    {}'.format(xopt1))
print ('Optimal function value:')
print ('    myfunc: {}'.format(fopt1))
#----------------------------------------------------------------------
##PySwarm
#print ('*'*65)
#print ('Example minimization of 4th-order banana function (with constraint)')
#def mycon(x):
#    x1 = x[0]
#    x2 = x[1]
#    return [-(x1 + 0.25)**2 + 0.75*x2]
#lb = min_ind
#ub = max_ind
#xopt2, fopt2 = pso(example_fitness, lb, ub, f_ieqcons=mycon, maxiter=200)
#print ('The optimum is at:')
#print ('    {}'.format(xopt2))
#print ('Optimal function value:')
#print ('    myfunc: {}'.format(fopt2))
#print ('    mycon : {}'.format(mycon(xopt2)))
#----------------------------------------------------------------------
##FST-PSO
#if (__name__ == "__main__"):
#  #dims = 22
#  FP = FuzzyPSO()
#  #FP.set_search_space( [[-3.0, 25.0]]*dims )
#  FP.set_search_space(min_max_ind)
#  FP.set_fitness(example_fitness)
#  #FSTPSO.set_swarm_size(100)
#  #bestpos, bestf = FSTPSO.solve_with_fstpso(max_iter=1000)
#  result =  FP.solve_with_fstpso(max_iter=10000)
#  print "------------------------"
#  print "Best solution:", result[0]
#  print "Whose fitness is:", result[1]
#----------------------------------------------------------------------

