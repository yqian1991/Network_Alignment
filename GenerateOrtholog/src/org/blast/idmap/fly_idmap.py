'''
Created on Apr 8, 2014

@author: yqian33
'''
def flyidmap(input, input2, output):
    
    count = 0
    #out = open('human_entrez2uniprot.txt','w')
    ins = open(input, 'r')
    out = open(output,'w')
    
    eachLine=list()
    for line in ins:
        if count>0 :
            eachLine = line.split('\t')
            item  = eachLine[5].split(':')
            #print eachLine[1]
            uniprotAC = findUniprotAC(input2, eachLine[1])
            #print uniprotAC
            if uniprotAC != 0:
                out.write(item[1]+'\t'+eachLine[1]+'\t'+uniprotAC+'\n')
        count = count + 1    
    ins.close()
    out.close()
    print 'Conversion Done'

def findUniprotAC(input, gid):  
    
    ins = open(input, 'r')
    
    str='0'
    flag=0
    for line in ins:
        eachline = line.split('\t')
        #print 'check:'+gid+'=='+eachline[1]
        #print gid==eachline[1]
        if gid == eachline[1] :
            #print 'hit:'+eachline[1]
            str = eachline[0]
            flag = 1
    #print count    
    ins.close()
    if flag == 1:
        #print 'find'
        return str
    else:
        return 0
        
def rawUniprotGid(input, output):
    ins = open(input, 'r')
    outs = open(output, 'w')
    for line in ins:
        eachline = line.split('\t')
        if eachline[2] != '':
            outs.write(eachline[0]+'\t'+eachline[2]+'\t\n')
    print 'raw data done'


def mapFlySeqName(input,input1, output):
    ins = open(input, 'r')
    ous = open(output, 'w')
    
    for line in ins:
        if line[0] == '>':
            tmpline = line.split('gene:')
            eachline = tmpline[1].split(' ')
            
            uniprotAC = substitiuteUniprotAC(input1, eachline[0])
            if uniprotAC != 0:
                ous.write('>'+uniprotAC)
        else:
            ous.write(line)
            
def substitiuteUniprotAC(input, fbid):  
    
    ins = open(input, 'r')
    
    str='0'
    flag=0
    for line in ins:
        eachline = line.split('\t')
        #print 'check:'+gid+'=='+eachline[1]
        #print gid==eachline[1]
        if fbid == eachline[0] :
            #print 'hit:'+eachline[1]
            str = eachline[2]
            flag = 1
    #print count    
    ins.close()
    if flag == 1:
        #print 'find'
        return str
    else:
        return 0
       
if __name__ == '__main__':
    
    infile1 = '/home/yqian33/APPI/dataset/ProteinSequences/Entrez_Drosophila_melanogaster.gene_info'
    infile2 = '/home/yqian33/APPI/dataset/ProteinSequences/DROME_7227_idmapping_selected.tab'
    outfile0 = '/home/yqian33/APPI/dataset/ProteinSequences/fly_gid_uniprot.txt'
    outfile = '/home/yqian33/APPI/dataset/ProteinSequences/flybase_uniprot.txt'
    #rawUniprotGid(infile2, outfile0)
    #flyidmap(infile1, outfile0, outfile)
    
    input = '/home/yqian33/APPI/dataset/ProteinSequences/Ensembl_Drosophila_melanogaster.BDGP5.75.pep.all.fa'
    #input = '/home/yqian33/APPI/dataset/ProteinSequences/fbgntest.fasta'
    output = '/home/yqian33/APPI/dataset/ProteinSequences/New_Ensembl_Drosophila_melanogaster.BDGP5.75.pep.all.fasta'
    
    import time
    start= time.clock()
    mapFlySeqName(input,outfile, output)
    end = time.clock()
    print end-start
    
    pass