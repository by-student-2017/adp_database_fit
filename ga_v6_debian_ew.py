import random
from deap import creator, base, tools, algorithms
import numpy
import numpy as np
import commands
import sys

print "Attention"
print "this code need Zhou04_ADP_EW"
print "  In the first time, need to command 'gfortran Zhou04_create_adp_ew.f -o Zhou04_ADP_EW'"
print "this code does not recommend element < Al"
print "all u parameters = 0.0"
print " "
#----------------------------------------------------------------------
file_tmp = 'ADP_code.tmp'
file_tmp_ew = 'ADP_code_ew.tmp'
file_inp = 'ADP_code'

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
commands.getoutput("chmod +x setinp")
commands.getoutput("./setinp")
commands.getoutput("chmod +x ./data2cfg/lmp_data2cfg")
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
print "read parameters from ADP_code.init"
nline = commands.getoutput("grep -n "+str(satom)+" ADP_code.init | head -1 | sed -e \"s/:.*//g\"")
print "read line: "+nline
check_satom = commands.getoutput("awk '{if(NR=="+str(nline)+"+0){print $1}}' ADP_code.init | head -1")
print "fit element: "+check_satom
# fitting parameters
x0  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+1){print $1}}' ADP_code.init | head -1"))
x1  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+2){print $1}}' ADP_code.init | head -1"))
x2  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+3){print $1}}' ADP_code.init | head -1"))
x3  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+4){print $1}}' ADP_code.init | head -1"))
x4  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+5){print $1}}' ADP_code.init | head -1"))
x5  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+6){print $1}}' ADP_code.init | head -1"))
x6  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+7){print $1}}' ADP_code.init | head -1"))
x7  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+8){print $1}}' ADP_code.init | head -1"))
x8  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+9){print $1}}' ADP_code.init | head -1"))
x9  = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+10){print $1}}' ADP_code.init | head -1"))
x10 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+11){print $1}}' ADP_code.init | head -1"))
x11 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+12){print $1}}' ADP_code.init | head -1"))
x12 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+13){print $1}}' ADP_code.init | head -1"))
x13 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+14){print $1}}' ADP_code.init | head -1"))
x14 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+15){print $1}}' ADP_code.init | head -1"))
x15 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+16){print $1}}' ADP_code.init | head -1"))
x16 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+17){print $1}}' ADP_code.init | head -1"))
x17 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+18){print $1}}' ADP_code.init | head -1"))
x18 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+19){print $1}}' ADP_code.init | head -1"))
x19 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+20){print $1}}' ADP_code.init | head -1"))
x20 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+23){print $1}}' ADP_code.init | head -1"))
x21 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+26){print $1}}' ADP_code.init | head -1"))
# u
x22 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+28){print $1}}' ADP_code.init | head -1"))
x23 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+29){print $1}}' ADP_code.init | head -1"))
x24 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+30){print $1}}' ADP_code.init | head -1"))
x25 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+31){print $1}}' ADP_code.init | head -1"))
x26 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+32){print $1}}' ADP_code.init | head -1"))
x27 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+33){print $1}}' ADP_code.init | head -1"))
x28 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+34){print $1}}' ADP_code.init | head -1"))
x29 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+35){print $1}}' ADP_code.init | head -1"))
x30 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+36){print $1}}' ADP_code.init | head -1"))
x31 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+37){print $1}}' ADP_code.init | head -1"))
x32 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+38){print $1}}' ADP_code.init | head -1"))
x33 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+39){print $1}}' ADP_code.init | head -1"))
x34 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+40){print $1}}' ADP_code.init | head -1"))
x35 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+41){print $1}}' ADP_code.init | head -1"))
x36 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+42){print $1}}' ADP_code.init | head -1"))
# w
x37 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+43){print $1}}' ADP_code.init | head -1"))
x38 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+44){print $1}}' ADP_code.init | head -1"))
x39 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+45){print $1}}' ADP_code.init | head -1"))
x40 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+46){print $1}}' ADP_code.init | head -1"))
x41 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+47){print $1}}' ADP_code.init | head -1"))
x42 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+48){print $1}}' ADP_code.init | head -1"))
x43 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+49){print $1}}' ADP_code.init | head -1"))
x44 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+50){print $1}}' ADP_code.init | head -1"))
x45 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+51){print $1}}' ADP_code.init | head -1"))
x46 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+52){print $1}}' ADP_code.init | head -1"))
x47 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+53){print $1}}' ADP_code.init | head -1"))
x48 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+54){print $1}}' ADP_code.init | head -1"))
x49 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+55){print $1}}' ADP_code.init | head -1"))
x50 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+56){print $1}}' ADP_code.init | head -1"))
x51 = float(commands.getoutput("awk '{if(NR=="+str(nline)+"+57){print $1}}' ADP_code.init | head -1"))

# u
fi0 = open(file_tmp,'r')
text0 = fi0.read().replace('urhoe',str(x22).replace("[","").replace("]",""))
text0 = text0.replace('un0',str(x23).replace("[","").replace("]",""))
text0 = text0.replace('un1',str(x24).replace("[","").replace("]",""))
text0 = text0.replace('un2',str(x25).replace("[","").replace("]",""))
text0 = text0.replace('un3',str(x26).replace("[","").replace("]",""))
text0 = text0.replace('u0',str(x27).replace("[","").replace("]",""))
text0 = text0.replace('u1',str(x28).replace("[","").replace("]",""))
text0 = text0.replace('u21',str(x29).replace("[","").replace("]",""))
text0 = text0.replace('u22',str(x30).replace("[","").replace("]",""))
text0 = text0.replace('u31',str(x31).replace("[","").replace("]",""))
text0 = text0.replace('u32',str(x32).replace("[","").replace("]",""))
text0 = text0.replace('ueta',str(x33).replace("[","").replace("]",""))
text0 = text0.replace('uep',str(x34).replace("[","").replace("]",""))
text0 = text0.replace('urhol',str(x35).replace("[","").replace("]",""))
text0 = text0.replace('urhoh',str(x36).replace("[","").replace("]",""))
with open(file_tmp_ew,'w') as f:
  print >> f, text0

x = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,
      x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51]
print "initial parameters: ", x

count = 0
#----------------------------------------------------------------------
creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

toolbox = base.Toolbox()

n_gene = 37 # number of parameters
min_ind = numpy.ones(n_gene) * -1.0
max_ind = numpy.ones(n_gene) *  1.0
for i in range(n_gene):
  #min_ind[i] = b1[i][0]
  #max_ind[i] = b1[i][1]
  min_ind[i] = float(x[i]) - float(x[i])*0.3
  max_ind[i] = float(x[i]) + float(x[i])*0.3
  print "search area of paramter "+str(i)+": "+str(min_ind[i])+" | "+str(max_ind[i])
#----------------------------------------------------------------------
def create_ind_uniform(min_ind, max_ind):
  ind = []
  for min, max in zip(min_ind, max_ind):
    ind.append(random.uniform(min, max))
  return ind
#----------------------------------------------------------------------
toolbox.register("create_ind", create_ind_uniform, min_ind, max_ind)
toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.create_ind)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
#----------------------------------------------------------------------
#def evalOneMax(individual):
#  return sum(individual),
#----------------------------------------------------------------------
def evalOneMax(individual):
  
  print "------------------------"
  global count
  count += 1
  print count

  fi = open(file_tmp_ew,'r')
  text = fi.read().replace('re',str(individual[0]).replace("[","").replace("]",""))
  text = text.replace('fe',str(individual[1]).replace("[","").replace("]",""))
  text = text.replace('Frhoe1',str(individual[2]).replace("[","").replace("]",""))
  text = text.replace('Frhoe2',str(individual[3]).replace("[","").replace("]",""))
  text = text.replace('alpha',str(individual[4]).replace("[","").replace("]",""))
  text = text.replace('beta',str(individual[5]).replace("[","").replace("]",""))
  text = text.replace('Ap',str(individual[6]).replace("[","").replace("]",""))
  text = text.replace('Bp',str(individual[7]).replace("[","").replace("]",""))
  text = text.replace('kappa',str(individual[8]).replace("[","").replace("]",""))
  text = text.replace('lambda',str(individual[9]).replace("[","").replace("]",""))
  text = text.replace('Fn0',str(individual[10]).replace("[","").replace("]",""))
  text = text.replace('Fn1',str(individual[11]).replace("[","").replace("]",""))
  text = text.replace('Fn2',str(individual[12]).replace("[","").replace("]",""))
  text = text.replace('Fn3',str(individual[13]).replace("[","").replace("]",""))
  text = text.replace('F0',str(individual[14]).replace("[","").replace("]",""))
  text = text.replace('F1',str(individual[15]).replace("[","").replace("]",""))
  text = text.replace('F2',str(individual[16]).replace("[","").replace("]",""))
  text = text.replace('F3',str(individual[17]).replace("[","").replace("]",""))
  text = text.replace('Feta',str(individual[18]).replace("[","").replace("]",""))
  text = text.replace('Fep',str(individual[19]).replace("[","").replace("]",""))
  text = text.replace('F4',str(individual[20]).replace("[","").replace("]",""))
  text = text.replace('Frhol',str(individual[21]).replace("[","").replace("]",""))
  # w
  text = text.replace('wrhoe',str(individual[22]).replace("[","").replace("]",""))
  text = text.replace('wn0',str(individual[23]).replace("[","").replace("]",""))
  text = text.replace('wn1',str(individual[24]).replace("[","").replace("]",""))
  text = text.replace('wn2',str(individual[25]).replace("[","").replace("]",""))
  text = text.replace('wn3',str(individual[26]).replace("[","").replace("]",""))
  text = text.replace('w0',str(individual[27]).replace("[","").replace("]",""))
  text = text.replace('w1',str(individual[28]).replace("[","").replace("]",""))
  text = text.replace('w21',str(individual[29]).replace("[","").replace("]",""))
  text = text.replace('w22',str(individual[30]).replace("[","").replace("]",""))
  text = text.replace('w31',str(individual[31]).replace("[","").replace("]",""))
  text = text.replace('w32',str(individual[32]).replace("[","").replace("]",""))
  text = text.replace('weta',str(individual[33]).replace("[","").replace("]",""))
  text = text.replace('wep',str(individual[34]).replace("[","").replace("]",""))
  text = text.replace('wrhol',str(individual[35]).replace("[","").replace("]",""))
  text = text.replace('wrhoh',str(individual[36]).replace("[","").replace("]",""))
  fi.close

  with open(file_inp,'w') as f:
    print >> f, text

  commands.getoutput("./Zhou04_ADP_1 < ADP.input")
  diffb  = commands.getoutput("cat diff.dat")
  if diffb == "nan" or abs(float(diffb)) >= 0.15/(1+float(count)/900):
    y = 999999.99999
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
      print "Lammps + PWscf calculation"
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
          fxl[i] = 99999999.99999
          fyl[i] = 99999999.99999
          fzl[i] = 99999999.99999
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
      if float(lmpe) == 0.0:
        lmpe = "99999999.99999"
  
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
      tdiffea = 99999999.99999
      tdiffp  = 99999999.99999
      tdifff  = 99999999.99999

  if error_flag1 != "" and error_flag2 == "" and error_flag3 == "":
    diffb  = commands.getoutput("cat diff.dat")
  else:
    diffb = 99999999.99999
  print "F boundary, diff: "+str(diffb)
  print "---------------"

  y = float(tdiffea)**2 + 1000*float(diffb)**2 + 0.0000002*abs(tdiffp)**2 + 0.0000010*abs(tdifff)**2

  print "Evaluate: ", y
  #print "Parameters: ", individual
  print "Parameters: x0 = "+"[ "+str(individual[0])+","+str(individual[1])+","+str(individual[2])+","+str(individual[3])+","+str(individual[4])+","+str(individual[5])+","+str(individual[6])+","+str(individual[7])+","+str(individual[8])+","+str(individual[9])+","+str(individual[10])+","+str(individual[11])+","+str(individual[12])+","+str(individual[13])+","+str(individual[14])+","+str(individual[15])+","+str(individual[16])+","+str(individual[17])+","+str(individual[18])+","+str(individual[19])+","+str(individual[20])+","+str(individual[21])+","+str(individual[22])+","+str(individual[23])+","+str(individual[24])+","+str(individual[25])+","+str(individual[26])+","+str(individual[27])+","+str(individual[28])+","+str(individual[29])+","+str(individual[30])+","+str(individual[31])+","+str(individual[32])+","+str(individual[33])+","+str(individual[34])+","+str(individual[35])+","+str(individual[36])+" ]"
  print "------------------------"
 
  return y,
#----------------------------------------------------------------------
def cxTwoPointCopy(ind1, ind2):
  size = len(ind1)
  cxpoint1 = random.randint(1, size)
  cxpoint2 = random.randint(1, size-1)
  if (cxpoint2 >= cxpoint1):
    cxpoint2 += 1
  else:
    cxpoint1, cxpoint2 = cxpoint2, cxpoint1

  ind1[cxpoint1:cxpoint2], ind2[cxpoint2:cxpoint2] = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()

  return ind1, ind2
#----------------------------------------------------------------------
def mutUniformDbl(individual, min_ind, max_ind, indpb):
  size = len(individual)
  for i, min, max in zip(xrange(size), min_ind, max_ind):
    if (random.random() < indpb):
      individual[i] = random.uniform(min, max)
  return indivisual,
#----------------------------------------------------------------------
toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
#----------------------------------------------------------------------
def main():
  random.seed(64)
  pop = toolbox.population(n=300)
  hof = tools.HallOfFame(1, similar=numpy.array_equal)
  stats = tools.Statistics(lambda ind: ind.fitness.values)
  stats.register("avg", numpy.mean)
  stats.register("std", numpy.std)
  stats.register("min", numpy.min)
  stats.register("max", numpy.max)
  algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=500, stats=stats, halloffame=hof)
  return pop, stats, hof
#----------------------------------------------------------------------
if (__name__ == "__main__"):
  main()
#----------------------------------------------------------------------

