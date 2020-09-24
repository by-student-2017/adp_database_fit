from scipy import optimize
import numpy
import numpy as np
import commands
import sys
#----------------------------------------------------------------------
file_tmp = 'ADP_code.tmp'
file_tmp_uw = 'ADP_code_uw.tmp'
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
if float(struct_list[3*ntemp+1]) <= 1073.0 :
  ntemp = ntemp + 1
  struct_list.append(1273.0)
  struct_list.append("L")
  struct_list.append(1.0)
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

fi0 = open(file_tmp,'r')
text0 = fi0.read().replace('re',str(x0).replace("[","").replace("]",""))
text0 = text0.replace('fe',str(x1).replace("[","").replace("]",""))
text0 = text0.replace('Frhoe1',str(x2).replace("[","").replace("]",""))
text0 = text0.replace('Frhoe2',str(x3).replace("[","").replace("]",""))
text0 = text0.replace('alpha',str(x4).replace("[","").replace("]",""))
text0 = text0.replace('beta',str(x5).replace("[","").replace("]",""))
text0 = text0.replace('Ap',str(x6).replace("[","").replace("]",""))
text0 = text0.replace('Bp',str(x7).replace("[","").replace("]",""))
text0 = text0.replace('kappa',str(x8).replace("[","").replace("]",""))
text0 = text0.replace('lambda',str(x9).replace("[","").replace("]",""))
text0 = text0.replace('Fn0',str(x10).replace("[","").replace("]",""))
text0 = text0.replace('Fn1',str(x11).replace("[","").replace("]",""))
text0 = text0.replace('Fn2',str(x12).replace("[","").replace("]",""))
text0 = text0.replace('Fn3',str(x13).replace("[","").replace("]",""))
text0 = text0.replace('F0',str(x14).replace("[","").replace("]",""))
text0 = text0.replace('F1',str(x15).replace("[","").replace("]",""))
text0 = text0.replace('F2',str(x16).replace("[","").replace("]",""))
text0 = text0.replace('F3',str(x17).replace("[","").replace("]",""))
text0 = text0.replace('Feta',str(x18).replace("[","").replace("]",""))
text0 = text0.replace('Fep',str(x19).replace("[","").replace("]",""))
text0 = text0.replace('F4',str(x20).replace("[","").replace("]",""))
text0 = text0.replace('Frhol',str(x21).replace("[","").replace("]",""))
with open(file_tmp_uw,'w') as f:
  print >> f, text0

z0 = [x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,
      x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51]
print "initial parameters: ",z0

count = 0
#----------------------------------------------------------------------
def f(z):
  
  print "------------------------"
  global count
  count += 1
  print count

  fi = open(file_tmp_uw,'r')
  # u
  text = fi.read().replace('urhoe',str(z[0]).replace("[","").replace("]",""))
  text = text.replace('un0',str(z[1]).replace("[","").replace("]",""))
  text = text.replace('un1',str(z[2]).replace("[","").replace("]",""))
  text = text.replace('un2',str(z[3]).replace("[","").replace("]",""))
  text = text.replace('un3',str(z[4]).replace("[","").replace("]",""))
  text = text.replace('u0',str(z[5]).replace("[","").replace("]",""))
  text = text.replace('u1',str(z[6]).replace("[","").replace("]",""))
  text = text.replace('u21',str(z[7]).replace("[","").replace("]",""))
  text = text.replace('u22',str(z[8]).replace("[","").replace("]",""))
  text = text.replace('u31',str(z[9]).replace("[","").replace("]",""))
  text = text.replace('u32',str(z[10]).replace("[","").replace("]",""))
  text = text.replace('ueta',str(z[11]).replace("[","").replace("]",""))
  text = text.replace('uep',str(z[12]).replace("[","").replace("]",""))
  text = text.replace('urhol',str(z[13]).replace("[","").replace("]",""))
  text = text.replace('urhoh',str(z[14]).replace("[","").replace("]",""))
  # w
  text = text.replace('wrhoe',str(z[15]).replace("[","").replace("]",""))
  text = text.replace('wn0',str(z[16]).replace("[","").replace("]",""))
  text = text.replace('wn1',str(z[17]).replace("[","").replace("]",""))
  text = text.replace('wn2',str(z[18]).replace("[","").replace("]",""))
  text = text.replace('wn3',str(z[19]).replace("[","").replace("]",""))
  text = text.replace('w0',str(z[20]).replace("[","").replace("]",""))
  text = text.replace('w1',str(z[21]).replace("[","").replace("]",""))
  text = text.replace('w21',str(z[22]).replace("[","").replace("]",""))
  text = text.replace('w22',str(z[23]).replace("[","").replace("]",""))
  text = text.replace('w31',str(z[24]).replace("[","").replace("]",""))
  text = text.replace('w32',str(z[25]).replace("[","").replace("]",""))
  text = text.replace('weta',str(z[26]).replace("[","").replace("]",""))
  text = text.replace('wep',str(z[27]).replace("[","").replace("]",""))
  text = text.replace('wrhol',str(z[28]).replace("[","").replace("]",""))
  text = text.replace('wrhoh',str(z[29]).replace("[","").replace("]",""))
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

    for itw in range(ntemp+1):
      if t == temp[itw]:
        wt = weig[i]

    tdiffea = tdiffea + float(diffea)*float(wt)
    tdiffp  = tdiffp  + float(diffp)*float(wt)
    tdifff  = tdifff  + float(difff)*float(wt)

  diffb  = commands.getoutput("cat diff.dat")
  print "F, u and w boundary, diff: "+str(diffb)
  print "---------------"

  y = float(tdiffea)**2 + 1000*float(diffb)**2 + 0.0000002*abs(tdiffp)**2 + 0.0000010*abs(tdifff)**2


  print "Evaluate: ", y
  print "Parameters: z0 = "+"[ "+str(z[0])+","+str(z[1])+","+str(z[2])+","+str(z[3])+","+str(z[4])+","+str(z[5])+","+str(z[6])+","+str(z[7])+","+str(z[8])+","+str(z[9])+","+str(z[10])+","+str(z[11])+","+str(z[12])+","+str(z[13])+","+str(z[14])+","+str(z[15])+","+str(z[16])+","+str(z[17])+","+str(z[18])+","+str(z[19])+","+str(z[20])+","+str(z[21])+","+str(z[22])+","+str(z[23])+","+str(z[24])+","+str(z[25])+","+str(z[26])+","+str(z[27])+","+str(z[28])+","+str(z[29])+" ]"
  print "------------------------"

  return y
#----------------------------------------------------------------------
res = optimize.minimize(f,z0,method='Nelder-Mead',options={'adaptive':True})
