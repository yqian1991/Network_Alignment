'''
Created on May 16, 2014

@author: yqian33
'''
def GOtermsByProtein(input, str):
    ins = open(input, 'r')
    
    print 'GO annotated to:'+str
    count=0
    list=[]
    for line in ins:
        if count>0:
            #print line
            item = line.split('\t')
            if item[1] == str:
                #print item[4]
                list.append(item[4])
        count+=1
    newlist = sorted(set(list)) 
    print   newlist
    print '--------------------------------------'
    return newlist 

def sharedTermsForProteins(list1, list2): 
    #list1=sorted(set(list1))
    #list2=sorted(set(list2))
    count = 0
    #print '---------------------------------------'
    print 'shared GO Terms:\n'
    for item1 in list1:
        for item2 in list2:
            if item1 == item2:
                print item1
                count +=1
    print 'count:'+str(count)           

def findGOtermsForlist():
    ins = open('/home/yqian33/APPI/software/AlignMCL-1.2/alignments/solutions/yeast.DIP-worm.DIP-yeast/99.txt', 'r')
    
    count=0
    total=0
    for line in ins:
        pname = line[:-1]
        total += 1
        glist = GOtermsByProtein(pname)
        if glist != []:
            count+=1
    
    print count
    print total
    
def mergeGOA(): 
    ins1=open('/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_fly', 'r')  
    ins2=open('/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_yeast', 'r')
    ous=open('/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_fly_yeast','w')
    
    count1=-1
    for line in ins1:
        ous.write(line)
        count1=count1+1
        
    count=0
    for line in ins2:
        if count>0:
            ous.write(line)
        count = count + 1
    print count1
    print count
    print count1+count
    print 'done'     
if __name__ == '__main__':
    #input1="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_yeast";
    #input2="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_fly"
    #list1 = GOtermsByProtein(input1, "P34761")
    #list2 = GOtermsByProtein(input2, "Q21559")
    #sharedTermsForProteins(list1, list2)
    #findGOtermsForlist()
    mergeGOA()
    pass