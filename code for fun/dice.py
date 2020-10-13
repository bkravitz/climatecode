import random

print "Enter dice you want to roll"
print "Examples:  1d6, 2d12, 2d8+1d4, etc."
print "You can also add modifiers, e.g., 3d6+5"

diceroll=raw_input("> ")

randarray=[]
mod=0
dicevals=0
results=[]
if '+' in diceroll:
    diceroll2=diceroll.split('+')
    for n in range(len(diceroll2)):
        temp=diceroll2[n]
        if 'd' in temp:
            temp2=temp.split('d')
            numtimes=int(temp2[0])
            dieval=int(temp2[1])
            while numtimes>0:
                randarray.append(dieval)
                numtimes=numtimes-1
        else:
            mod=mod+int(temp)
else:
    if 'd' in diceroll:
        diceroll2=diceroll.split('d')
        numtimes=int(diceroll2[0])
        dieval=int(diceroll2[1])
        while numtimes>0:
            randarray.append(dieval)
            numtimes=numtimes-1
    else:
        mod=mod+int(diceroll)

for k in range(len(randarray)):
    rollval=random.randrange(1,randarray[k]+1,1)
    if randarray[k]>9:
        spaces='  '
    else:
        spaces='   '
    outputstr='1d'+str(randarray[k])+':'+spaces+str(rollval)
    results.append(outputstr)
    dicevals=dicevals+rollval

totalvals=dicevals+mod

print "Rolls:"
for j in range(len(results)):
    print results[j]
print "\nTotal: " + str(totalvals)