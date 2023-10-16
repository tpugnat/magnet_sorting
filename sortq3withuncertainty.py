from numpy import *
from numpy.random import normal
from random import *
#import matplotlib.pyplot as plt

def sortOnSum(val): 
    return sum(val) 



Nq=8     # Number of individual Q2, 2 ips, 2 quads per ip side
FullRangeT=50  # Full range of integrated gradient error in units
measerror = 2  # Random error of measurement

Alldiffs=[]
Allaves=[]
Allgroups=[]

for r in range(1):
 systerror = 10*2*(random()-0.5) # systematic error of meas
 print "systerror", systerror
 q=[]
 m=[]
 for i in range(Nq):
    q.append(FullRangeT*(random()-0.5)) # generate real errors
    m.append(measerror*2*(random()-0.5)) # generate real errors
    

 qm=zip(q,m)
 #print qm
 #qm.sort(key=sortOnSum) # Sort according to measurement
 #print qm
 print "Actual magnet errors of paired magnets"
 for i in range(0,len(qm)-1,2):
  print qm[i][0], qm[i+1][0]

 print "Believed magnet errors of paired magnets"
 for i in range(0,len(qm)-1,2):
  print sum(qm[i]) + systerror, sum(qm[i+1]) + systerror

 print "Pairs powered as"
 for i in range(0,len(qm)-1,2):
  print (sum(qm[i]) + systerror + sum(qm[i+1]) + systerror)/2

 print "Errors to put in MADX"
 error_file = open('errors_Q3.madx','w')
 #LHC
 #magnet_list = ["MQXA.3L1", "MQXA.3R1", "MQXA.3L5", "MQXA.3R5"]
 #HL-LHC
 magnet_list = ["MQXFA.A3L1", "MQXFA.B3L1", "MQXFA.A3R1", "MQXFA.B3R1", 
 "MQXFA.A3L5", "MQXFA.B3L5", "MQXFA.A3R5", "MQXFA.B3R5"]
 qgrouped=[]
 for i in range(0,len(qm)-1,2):
  pairpower=(sum(qm[i]) + systerror + sum(qm[i+1]) + systerror)/2
  qgrouped.append( [-qm[i][0] + pairpower, -qm[i+1][0] + pairpower])
  print qgrouped[-1]
  #print>>error_file, qgrouped[-1]
  for j in range(2):
	  print>>error_file, 'select, flag=error, clear;'
	  print>>error_file, 'select, flag=error, PATTERN=', magnet_list[i+j],';'
	  print>>error_file, 'Efcomp, radius = Rr, order= 1, dknr:={0,1E-4*(',qgrouped[-1][j],')};'
	  #print>>error_file, 'B2r =',qgrouped[-1][j],';'
	  #print>>error_file, 'exec SetEfcomp_Q;'
	  print>>error_file, ''

 
 #print qgrouped
 diff=[]
 for pair in qgrouped:
    diff.append(abs(pair[1]-pair[0]))

 Allaves.append(average(qgrouped))
 Alldiffs = Alldiffs + list(diff)
 Allgroups.append(array(qgrouped)-Allaves[-1])


Allq=array(Allgroups).flatten()

print std(Allq),average(Alldiffs), std(Alldiffs), max(Alldiffs), average(Allaves)

#plt.hist(Allq, cumulative=False)
#plt.show()
#plt.hist(Alldiffs,cumulative=False)
#plt.show()
#plt.hist(Alldiffs,cumulative=True)
#plt.show()
#plt.hist(Allaves,cumulative=False)
#plt.show()

# No pop
# 0.272840491993 0.24857240342 2.21996240797
#POP
#0.185596047974 0.150069713415 1.03630474542
