from matplotlib import pyplot as plt
import numpy as np
import argparse

#loading data file and parsing line by line creating a list of mutations and a list of positions 
spl =[] #mutation list 
pos=[] #position list 
#n=5
#rep=100

parser = argparse.ArgumentParser(
                    prog = 'discoal_data.py',
                    description = 'plot hardsweep data to compare with msprime',
                    epilog = '-txt')
parser.add_argument('-txt', type=str, default="discoal_default.txt")
args=parser.parse_args()
file=args.txt

with open(file) as f:
    L=f.readline()
    param=L.split()
    n=int(param[1])
    rep=int(param[2])
    print(n)
    
    L=f.readline()
    L=f.readline()
    prev=""
    reps=0 #initially starting at 0 reps
    pos=[]
    while reps<=rep+1:
        spl_sub=[]
        cnt=0
        while cnt<n+2:
            L=f.readline()
            if "//" in L:
                break
            elif "position" in L:
                pos_list=(L[L.index(":")+2:-2].split(" "))
                pos_list=[float(x) for x in pos_list[:-1]]
                pos.append(pos_list)

                
            elif "discoal" not in L and "segsites" not in L:
                #print(L)
                add=list(L)
                res = [eval(j) for j in add[:-1]] #converting form str to int
                spl_sub.append(res)
                #print(spl_sub)
                #L=f.readline()  
                #print(L)
            cnt+=1
        reps+=1 #incrementing to keep track of reps
        #print(reps)
        spl.append(spl_sub)
del spl[0::2]
#print(spl)
#print(pos)

#creating an array of reduced rows for all reps with an array per rep
sumarr=[]
for i in range(len(spl)):
    sumarr.append([sum(x) for x in zip(*spl[i])])
#print(sumarr)

#getting array of heterozygosity for each rep
total_het=[]
for i in range (len(sumarr)):
    
    #getting a list of heterozygosity at each position
    heterozygosity=[]
    for j in range(len(sumarr[i])):
        val=sumarr[i][j]*2*(n-sumarr[i][j])/(n*(n-1))
        heterozygosity.append(val)
    total_het.append(heterozygosity)
#print(total_het)
#print(type(pos[0][0]))



#creating a list of windows 
win=[0]
win_het=[]
for i in range(1,12):
    win.append(i/11)
#print(win)



for i in range(len(pos)):
    win_het.append(list(0 for i in range(len(win)-1)))
#print(win_het)



#print(total_het)
for j in range(len(pos)):
    for i in range(len(pos[j])):
        for k in range(len(win)-1):
            if(pos[j][i]>=win[k] and pos[j][i]<win[k+1]):
                win_het[j][k]+=total_het[j][i]/((10**6)/11)
#print(win_het)
x=[i for i in range(1,12)]

win_avg=np.average(win_het,axis=0)

#print(win_avg)
#plt.plot(x,win_avg, marker='o')

win_box=[]
for i in range(len(win_het[0])):
    temp=[]
    for j in range(len(win_het)):
        temp.append(win_het[j][i])
    win_box.append(temp)
                   
positions=range(1,len(win_box)+1)               
                   
#plt.boxplot(win_box, positions)
#plt.errorbar(x, win_avg, linestyle='None', marker='^')

with open('discoal_data.npy', 'wb') as f:
    np.save(f, win_avg, allow_pickle=True)


