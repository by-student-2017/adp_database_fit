from scipy import optimize
import numpy
import numpy as np
import commands
import sys
#----------------------------------------------------------------------
file_tmp = 'ADP_code.tmp'
file_inp = 'ADP_code'

satom = commands.getoutput("grep \"atomtype\" ADP.input | sed -e \"s/.*=//\" -e \"s/'//g\"")

commands.getoutput("chmod +x setinp")
commands.getoutput("./setinp")
#----------------------------------------------------------------------

x0 = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,
      x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,
      x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51]
print "initial parameters: ",x0

#----------------------------------------------------------------------
def f(x):
  
  print "------------------------"

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

  y = 0.0

  print "made EAM potential"
  print "  "+str(satom)+"_Zhou04.adp"
  print "------------------------"

  return y
#----------------------------------------------------------------------
res = f(x0)
#res = optimize.minimize(f,x0,method='Nelder-Mead',options={'adaptive':True})
#res = optimize.minimize(f,x0,method='Nelder-Mead')
#res = optimize.minimize(f,x0,method='TNC')
#res = optimize.minimize(f,x0,method='Powell')
#res = optimize.minimize(f,x0,method='BFGS')
