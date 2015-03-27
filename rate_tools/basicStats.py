from RptMatrix import *
import operator
import pickle
from asymm_tools import compute_d
from asymm_tools import compute_P

rpt_list = load_psm("/Users/mcallimb/cache/human/hg18/seq/rmsk/chr9.psm")
counts = {}
mats = {}
f = open('./workfile.txt', 'w')
for i in rpt_list:
    if i.class_name in counts:
        counts[i.class_name]+=1
        mats[i.class_name] += i.M
        f.write(i.class_name +": " + str(i.start) + ";" + str(i.finish) +"\n")
    else:
        counts[i.class_name] = 1
        mats[i.class_name] = i.M
sorted_x = sorted(counts.items(), key=operator.itemgetter(1),reverse=True)
sum = 0
for x in sorted_x:
        print x
        print mats[x[0]]
        print compute_d(compute_P(mats[x[0]]))
        if x[0][:3] == "Alu":
            sum += x[1] 
print sum