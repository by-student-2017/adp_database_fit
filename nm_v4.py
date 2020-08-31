from scipy import optimize
import numpy
import numpy as np
import commands
import sys
#----------------------------------------------------------------------
file_tmp = 'ADP_code.tmp'
file_inp = 'ADP_code'

cif2cell_adress = "cif2cell"

commands.getoutput("setenv OMP_NUM_THREADS 1")
num_core = commands.getoutput("grep 'core id' /proc/cpuinfo | sort -u | wc -l")
#lammps_adress = "mpirun -np "+str(num_core)+" --allow-run-as-root lmp"
#pwscf_adress = "mpirun -np "+str(num_core)+" --allow-run-as-root pw.x"
#lammps_adress = "mpirun -np "+str(num_core)+" lmp"
#pwscf_adress = "mpirun -np "+str(num_core)+" pw.x"
lammps_adress = "mpirun -np 2 lmp"
pwscf_adress = "mpirun -np 2 pw.x"

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

#temp = [77,423,1073]
temp = [77]
for t in temp:
  commands.getoutput("cp in.lmp in.lmp_"+str(t)+"K")
  commands.getoutput("sed -i 's/YYYY/"+str(t)+"/' in.lmp_"+str(t)+"K")
  commands.getoutput("cp data.in data.in_"+str(t)+"K")
# if you would set other structures, e.g.,
commands.getoutput("cp ./data/data.in.DC data.in_423K")
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
x0 = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,
      x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,
      x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51]
print "initial parameters: ",x0

count = 0
#----------------------------------------------------------------------
def f(x):
  
  print "------------------------"
  global count
  count += 1
  print count

  fi = open(file_tmp,'r')
  text = fi.read().replace('re',str(x[0]).replace("[","").replace("]",""))
  text = text.replace('fe',str(x[1]).replace("[","").replace("]",""))
  text = text.replace('Frhoe1',str(x[2]).replace("[","").replace("]",""))
  text = text.replace('Frhoe2',str(x[3]).replace("[","").replace("]",""))
  text = text.replace('alpha',str(x[4]).replace("[","").replace("]",""))
  text = text.replace('beta',str(x[5]).replace("[","").replace("]",""))
  text = text.replace('Ap',str(x[6]).replace("[","").replace("]",""))
  text = text.replace('Bp',str(x[7]).replace("[","").replace("]",""))
  text = text.replace('kappa',str(x[8]).replace("[","").replace("]",""))
  text = text.replace('lambda',str(x[9]).replace("[","").replace("]",""))
  text = text.replace('Fn0',str(x[10]).replace("[","").replace("]",""))
  text = text.replace('Fn1',str(x[11]).replace("[","").replace("]",""))
  text = text.replace('Fn2',str(x[12]).replace("[","").replace("]",""))
  text = text.replace('Fn3',str(x[13]).replace("[","").replace("]",""))
  text = text.replace('F0',str(x[14]).replace("[","").replace("]",""))
  text = text.replace('F1',str(x[15]).replace("[","").replace("]",""))
  text = text.replace('F2',str(x[16]).replace("[","").replace("]",""))
  text = text.replace('F3',str(x[17]).replace("[","").replace("]",""))
  text = text.replace('Feta',str(x[18]).replace("[","").replace("]",""))
  text = text.replace('Fep',str(x[19]).replace("[","").replace("]",""))
  text = text.replace('F4',str(x[20]).replace("[","").replace("]",""))
  text = text.replace('Frhol',str(x[21]).replace("[","").replace("]",""))
  # u
  text = text.replace('urhoe',str(x[22]).replace("[","").replace("]",""))
  text = text.replace('un0',str(x[23]).replace("[","").replace("]",""))
  text = text.replace('un1',str(x[24]).replace("[","").replace("]",""))
  text = text.replace('un2',str(x[25]).replace("[","").replace("]",""))
  text = text.replace('un3',str(x[26]).replace("[","").replace("]",""))
  text = text.replace('u0',str(x[27]).replace("[","").replace("]",""))
  text = text.replace('u1',str(x[28]).replace("[","").replace("]",""))
  text = text.replace('u21',str(x[29]).replace("[","").replace("]",""))
  text = text.replace('u22',str(x[30]).replace("[","").replace("]",""))
  text = text.replace('u31',str(x[31]).replace("[","").replace("]",""))
  text = text.replace('u32',str(x[32]).replace("[","").replace("]",""))
  text = text.replace('ueta',str(x[33]).replace("[","").replace("]",""))
  text = text.replace('uep',str(x[34]).replace("[","").replace("]",""))
  text = text.replace('urhol',str(x[35]).replace("[","").replace("]",""))
  text = text.replace('urhoh',str(x[36]).replace("[","").replace("]",""))
  # w
  text = text.replace('wrhoe',str(x[37]).replace("[","").replace("]",""))
  text = text.replace('wn0',str(x[38]).replace("[","").replace("]",""))
  text = text.replace('wn1',str(x[39]).replace("[","").replace("]",""))
  text = text.replace('wn2',str(x[40]).replace("[","").replace("]",""))
  text = text.replace('wn3',str(x[41]).replace("[","").replace("]",""))
  text = text.replace('w0',str(x[42]).replace("[","").replace("]",""))
  text = text.replace('w1',str(x[43]).replace("[","").replace("]",""))
  text = text.replace('w21',str(x[44]).replace("[","").replace("]",""))
  text = text.replace('w22',str(x[45]).replace("[","").replace("]",""))
  text = text.replace('w31',str(x[46]).replace("[","").replace("]",""))
  text = text.replace('w32',str(x[47]).replace("[","").replace("]",""))
  text = text.replace('weta',str(x[48]).replace("[","").replace("]",""))
  text = text.replace('wep',str(x[49]).replace("[","").replace("]",""))
  text = text.replace('wrhol',str(x[50]).replace("[","").replace("]",""))
  text = text.replace('wrhoh',str(x[51]).replace("[","").replace("]",""))
  fi.close

  with open(file_inp,'w') as f:
    print >> f, text

  commands.getoutput("./Zhou04_ADP_1 < ADP.input")

  tdiffea = 0.0
  tdiffp  = 0.0
  tdifff  = 0.0
  for t in temp:
    print "---------------"
    print "Temperature: "+str(t)+" [K]"
    if count > 9000 or count % int(3000*2.718**(-count/3000)+1) == 1:
      commands.getoutput("mv data.in_"+str(t)+"K data.in")
      natom = commands.getoutput("awk '{if($2==\"atoms\"){print $1}}' data.in")
      commands.getoutput(lammps_adress+" < in.lmp_"+str(t)+"K")
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

    # stress = pressure
    pxxl = commands.getoutput("awk '{if($1==\"pxxl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    pyyl = commands.getoutput("awk '{if($1==\"pyyl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    pzzl = commands.getoutput("awk '{if($1==\"pzzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    pxyl = commands.getoutput("awk '{if($1==\"pxyl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    pxzl = commands.getoutput("awk '{if($1==\"pxzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    pyzl = commands.getoutput("awk '{if($1==\"pyzl\"){printf \"%10.8f\",$3*7.4028083e-11}}' log.lammps")
    pxxp = commands.getoutput("awk '{if($1==\"#S\"){print $2}}' config_"+str(t)+"K")
    pyyp = commands.getoutput("awk '{if($1==\"#S\"){print $3}}' config_"+str(t)+"K")
    pzzp = commands.getoutput("awk '{if($1==\"#S\"){print $4}}' config_"+str(t)+"K")
    pxyp = commands.getoutput("awk '{if($1==\"#S\"){print $5}}' config_"+str(t)+"K")
    pxzp = commands.getoutput("awk '{if($1==\"#S\"){print $6}}' config_"+str(t)+"K")
    pyzp = commands.getoutput("awk '{if($1==\"#S\"){print $7}}' config_"+str(t)+"K")
    diffpxx = (float(pxxl) - float(pxxp))/(float(pxxp)+0.000000101)*100.0/6.0
    diffpyy = (float(pyyl) - float(pyyp))/(float(pyyp)+0.000000101)*100.0/6.0
    diffpzz = (float(pzzl) - float(pzzp))/(float(pzzp)+0.000000101)*100.0/6.0
    diffpxy = (float(pxyl) - float(pxyp))/(float(pxyp)+0.000000101)*100.0/6.0
    diffpxz = (float(pxzl) - float(pxzp))/(float(pxzp)+0.000000101)*100.0/6.0
    diffpyz = (float(pyzl) - float(pyzp))/(float(pyzp)+0.000000101)*100.0/6.0
    diffp = abs(diffpxx) + abs(diffpyy) + abs(diffpzz) + abs(diffpxy) + abs(diffpxz) + abs(diffpyz)
    print "lammps: "+str(pxxl)+", "+str(pyyl)+", "+str(pzzl)+", "+str(pxyl)+", "+str(pxzl)+", "+str(pyzl)+" [eV/A^3]"
    print "PWscf:  "+str(pxxp)+", "+str(pyyp)+", "+str(pzzp)+", "+str(pxyp)+", "+str(pxzp)+", "+str(pyzp)+" [eV/A^3]"
    print "P diff (%): "+str(diffp)
    print "---------------"

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

    tdiffea = tdiffea + float(diffea)
    tdiffp  = tdiffp  + float(diffp)
    tdifff  = tdifff  + float(difff)

  rhoin  = float(x[2])*float(x[21])
  rhoout = float(x[2])*1.15
  print "---------------"
  print "F boundary 1, rho: "+str(rhoin)
  print "F boundary 2, rho: "+str(x[2])
  print "F boundary 3, rho: "+str(rhoout)
  commands.getoutput("cp "+satom+"_Zhou04.adp"+" Xx_Zhou04.adp")
  commands.getoutput("./plot_adp")
  rhoin1  = commands.getoutput("cat F.plt | awk '{if($1<"+str(rhoin)+"){print $2}}' | tail -2 | head -1")
  rhoin2  = commands.getoutput("cat F.plt | awk '{if($1>"+str(rhoin)+"){print $2}}' | head -2 | tail -1")
  rhoe1   = commands.getoutput("cat F.plt | awk '{if($1<"+str(x[2])+"){print $2}}' | tail -2 | head -1")
  rhoe2   = commands.getoutput("cat F.plt | awk '{if($1>"+str(x[2])+"){print $2}}' | head -2 | tail -1")
  rhoout1 = commands.getoutput("cat F.plt | awk '{if($1<"+str(rhoout)+"){print $2}}' | tail -2 | head -1")
  rhoout2 = commands.getoutput("cat F.plt | awk '{if($1>"+str(rhoout)+"){print $2}}' | head -2 | tail -1")
  print "near F boundary 1, F: "+str(rhoin1)+" | "+str(rhoin2)+" | diff "+str(float(rhoin1) - float(rhoin2))
  print "near F boundary 2, F: "+str(rhoe1)+" | "+str(rhoe2)+" | diff "+str(float(rhoe1) - float(rhoe2))
  print "near F boundary 3, F: "+str(rhoout1)+" | "+str(rhoout2)+" | diff "+str(float(rhoout1) - float(rhoout2))
  print "---------------"

  urhoin  = float(x[22])*float(x[35])
  urhoout = float(x[22])*float(x[36])
  print "---------------"
  print "u boundary 1, rho: "+str(urhoin)
  print "u boundary 2, rho: "+str(x[22])
  print "u boundary 3, rho: "+str(urhoout)
  urhoin1  = commands.getoutput("cat u.plt | awk '{if($1<"+str(urhoin)+"){print $2}}' | tail -2 | head -1")
  urhoin2  = commands.getoutput("cat u.plt | awk '{if($1>"+str(urhoin)+"){print $2}}' | head -2 | tail -1")
  urhoe1   = commands.getoutput("cat u.plt | awk '{if($1<"+str(x[22])+"){print $2}}' | tail -2 | head -1")
  urhoe2   = commands.getoutput("cat u.plt | awk '{if($1>"+str(x[22])+"){print $2}}' | head -2 | tail -1")
  urhoout1 = commands.getoutput("cat u.plt | awk '{if($1<"+str(urhoout)+"){print $2}}' | tail -2 | head -1")
  urhoout2 = commands.getoutput("cat u.plt | awk '{if($1>"+str(urhoout)+"){print $2}}' | head -2 | tail -1")
  print "near u boundary 1, u: "+str(urhoin1)+" | "+str(urhoin2)+" | diff "+str(float(urhoin1) - float(urhoin2))
  print "near u boundary 2, u: "+str(urhoe1)+" | "+str(urhoe2)+" | diff "+str(float(urhoe1) - float(urhoe2))
  print "near u boundary 3, u: "+str(urhoout1)+" | "+str(urhoout2)+" | diff "+str(float(urhoout1) - float(urhoout2))
  print "---------------"

  wrhoin  = float(x[37])*float(x[50])
  wrhoout = float(x[37])*float(x[51])
  print "---------------"
  print "w boundary 1, rho: "+str(wrhoin)
  print "w boundary 2, rho: "+str(x[37])
  print "w boundary 3, rho: "+str(wrhoout)
  wrhoin1  = commands.getoutput("cat w.plt | awk '{if($1<"+str(wrhoin)+"){print $2}}' | tail -2 | head -1")
  wrhoin2  = commands.getoutput("cat w.plt | awk '{if($1>"+str(wrhoin)+"){print $2}}' | head -2 | tail -1")
  wrhoe1   = commands.getoutput("cat w.plt | awk '{if($1<"+str(x[37])+"){print $2}}' | tail -2 | head -1")
  wrhoe2   = commands.getoutput("cat w.plt | awk '{if($1>"+str(x[37])+"){print $2}}' | head -2 | tail -1")
  wrhoout1 = commands.getoutput("cat w.plt | awk '{if($1<"+str(wrhoout)+"){print $2}}' | tail -2 | head -1")
  wrhoout2 = commands.getoutput("cat w.plt | awk '{if($1>"+str(wrhoout)+"){print $2}}' | head -2 | tail -1")
  print "near w boundary 1, w: "+str(wrhoin1)+" | "+str(wrhoin2)+" | diff "+str(float(wrhoin1) - float(wrhoin2))
  print "near w boundary 2, w: "+str(wrhoe1)+" | "+str(wrhoe2)+" | diff "+str(float(wrhoe1) - float(wrhoe2))
  print "near w boundary 3, w: "+str(wrhoout1)+" | "+str(wrhoout2)+" | diff "+str(float(wrhoout1) - float(wrhoout2))
  print "---------------"

  y = abs(tdiffea)**2 + 1000*abs(float(rhoin1) - float(rhoin2))**2 + 1000*abs(float(rhoe1) - float(rhoe2))**2 + 1000*abs(float(rhoout1) - float(rhoout2))**2 + 0.0000002*abs(tdiffp)**2 + 0.0000010*abs(tdifff)**2  + 1000*abs(float(urhoin1) - float(urhoin2))**2 + 1000*abs(float(urhoe1) - float(urhoe2))**2 + 1000*abs(float(urhoout1) - float(urhoout2))**2 + 1000*abs(float(wrhoin1) - float(wrhoin2))**2 + 1000*abs(float(wrhoe1) - float(wrhoe2))**2 + 1000*abs(float(wrhoout1) - float(wrhoout2))**2

  print "Evaluate: ", y
  #print "Parameters: ", x
  print "Parameters: x0 = "+"[ "+str(x[0])+","+str(x[1])+","+str(x[2])+","+str(x[3])+","+str(x[4])+","+str(x[5])+","+str(x[6])+","+str(x[7])+","+str(x[8])+","+str(x[9])+","+str(x[10])+","+str(x[11])+","+str(x[12])+","+str(x[13])+","+str(x[14])+","+str(x[15])+","+str(x[16])+","+str(x[17])+","+str(x[18])+","+str(x[19])+","+str(x[20])+","+str(x[21])+","+str(x[22])+","+str(x[23])+","+str(x[24])+","+str(x[25])+","+str(x[26])+","+str(x[27])+","+str(x[28])+","+str(x[29])+","+str(x[30])+","+str(x[31])+","+str(x[32])+","+str(x[33])+","+str(x[34])+","+str(x[35])+","+str(x[36])+","+str(x[37])+","+str(x[38])+","+str(x[39])+","+str(x[40])+","+str(x[41])+","+str(x[42])+","+str(x[43])+","+str(x[44])+","+str(x[45])+","+str(x[46])+","+str(x[47])+","+str(x[48])+","+str(x[49])+","+str(x[50])+","+str(x[51])+" ]"
  print "------------------------"

  return y
#----------------------------------------------------------------------
res = optimize.minimize(f,x0,method='Nelder-Mead',options={'adaptive':True})
#res = optimize.minimize(f,x0,method='Nelder-Mead')
#res = optimize.minimize(f,x0,method='TNC')
#res = optimize.minimize(f,x0,method='Powell')
#res = optimize.minimize(f,x0,method='BFGS')
