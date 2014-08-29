'''
Created on Mar 18, 2014

@author: yqian33
'''
 
 
    
'''
Generate   Ortholog file For Two Species  
'''

#########################################################################
#                    Generate   Ortholog file For Two Species           #
#########################################################################
def extractIDFromDIPOT(input, output, output1):
    
    import csv
    
    #ous= open('/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/fly_worm/fly_worm_ort.txt','w') 
    ous= open(output1, 'w')
    #ous1 = open('/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/fly_worm/worm_geneid.txt', 'w')
    ous1 = open(output, 'w')
    ous2 = open('/home/yqian33/workspace/PPIDataProcess/GenerateOrtholog/src/org/blast/evalue/yeast_human/temp/human_mgi.txt', 'w')
    
    count = 0;
    #/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/fly_worm/fly_worm_ort.xls
    with open(input, 'r') as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            if count > 0:
                if line[7] != '0':
                    ous.write(line[0]+'\t'+line[4]+'\t'+line[5]+'\t'+str( round((float(line[8])/float(line[7])), 2) )+'\n')
                    ous1.write(line[4]+'\n')
                    ous2.write(line[5]+'\n')
            count = count + 1
    ous.close()
    ous1.close()
    ous2.close()
    print 'done'

def proteinMask(input, output):
    ins = open(input, 'r')
    ous = open(output, 'w')
    
    plist=[]
    count = 0
    for line in ins:
        flag = idExisted(plist, line[:-1])
        if flag == 0:
            plist.append(line[:-1])
            ous.write(line)
            count += 1
    ins.close()
    ous.close()
    
def idExisted(plist, str):
    flag = 0
    
    num = len(plist)
    for count in range(0, num):
        if str == plist[count]:
            flag = 1
            return flag
        
    return flag

def mapGid2Uac(input, mapinput, output):
    
    #ins = open( '/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/human_fly/fly_geneid.txt','r')
    ins = open(input,'r')
    
    #ous = open('/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/human_fly/fly_gid_uac.txt', 'w')
    ous = open(output, 'w')
    count = 0
    
    for line in ins:
        
        if line[-1] == '\n':
            #print line[:-1]
            #str = searchMapFile('/home/yqian33/APPI/ProSeqBla/idmapfile/DROME_7227_idmapping_selected.tab', line[:-1])
            str = searchMapFile(mapinput, line[:-1])
            #print str
        else:
            #print line
            str = searchMapFile(mapinput, line)
            #print str
        
        if str != 0:
            stritem = str.split('\t')
            for item in stritem[:-1]:
                ous.write(line[:-1]+'\t'+item+'\n')
                count = count + 1
    
    print count
    ins.close()
    ous.close()
    print 'done'

def searchMapFile(input, str):
    
    ins = open(input, 'r')
    
    result =''
    print str
    for line in ins:
        
        item = line.split('\t')
        #print 'line'+item[2]
        if ';' in item[2]:
            iitem = item[2].split(';')
            for col in iitem:
                #print col
                if col == str:
                    #print col
                    result += item[0]+'\t'
        else:
            if item[2] == str:
                #print item[2]
                #print item[2]
                result += item[0] +'\t'
            #return item[0]
    
    ins.close()
    if result != '':
        #print 'good'+result
        return result
    else:
        return 0

class Protein:
    
    def __init__(self, entrezid='', uniprotac=''):
        self.entrezid = entrezid
        self.uniprotac = uniprotac
    
    def __del__(self):
        pass
    
def repeatMask(input, output):
    
    #ins = open('/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/human_fly/fly_gid_uac_clean.txt', 'r')
    #ous = open("/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/human_fly/fly_gid_uac_test.txt", 'w')
    ins=open(input, 'r')
    ous=open(output, 'w')
    plist=[]
    flag = 0
    count = 0
    for line in ins:
        str = line.split('\t')
        p=Protein(str[0], str[1][:-1])
        flag = existed(plist, p)
        if flag == 0:
            plist.append(p)
            ous.write(str[0]+'\t'+str[1][:-1]+'\n')
            count +=1
    ins.close()
    ous.close()
    print count
    print 'done'
            
def existed(plist, p):   
    flag = 0
    num = len(plist)
    for count in range(0, num):
        po = plist[count]
        if po.entrezid == p.entrezid and po.uniprotac == p.uniprotac:
            flag = 1
            return flag
    
    return flag          
    
def map(input, output, midinput):
    
    #ins = open('/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/fly_worm/fly_worm_ort.txt','r') 
    ins = open(input, 'r')
    #ous = open('/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/fly_worm/fly_worm.ort','w')
    ous = open(output, 'w')
    
    count = 0
    for line in ins:
        item = line.split('\t')
        #item[0], item[1]=geneid, item[2]=sgdid, item[3]=score
        geneidAC = searchAC(midinput, item[1])
        #sgdidAC = searchAC('/home/yqian33/workspace/PPIDataProcess/BLASTPforEvalue/src/org/blast/evalue/fly_worm/worm_wbid_uniprotAC.txt', item[2])
        
        if geneidAC != 0:
            ous.write(item[0]+'\t'+geneidAC[:-1]+'\t'+item[3])
            count = count + 1
            #print item[3]
        #elif sgdidAC != 0:
        #   count = count + 1
            #ous.write(item[0]+'\t'+sgdidAC[:-1]+'\t'+item[3])
        else:
            print item[1]+' '+item[2]+'***No hit***'
    
    ins.close()
    ous.close()
    print count
    print 'done'

def searchAC(input, str):
    
    ins = open(input, 'r')
    
    flag = 0
    for line in ins:
        item = line.split('\t')
        if item[0] == str:
            flag = 1
            return item[1]
        else:
            flag = 0
    
    ins.close()
    if flag == 0:
        return 0 

def idmapping():
    input='/home/yqian33/workspace/PPIDataProcess/GenerateOrtholog/src/org/blast/evalue/yeast_human/yeast_human_ort.xls'
    output='/home/yqian33/workspace/PPIDataProcess/GenerateOrtholog/src/org/blast/evalue/yeast_human/temp/human_geneid.txt'
    
    #ort file after first clean(id is uniprot ac and gene name)
    output1='/home/yqian33/workspace/PPIDataProcess/GenerateOrtholog/src/org/blast/evalue/yeast_human/temp/yeast_human_ort.txt'
    mapinput ='/home/yqian33/APPI/ProSeqBla/idmapfile/HUMAN_9606_idmapping_selected.tab'
    
    print 'Extract protein gene id for species2...'
    extractIDFromDIPOT(input, output, output1)
    print 'done'
    
    maskidout = '/home/yqian33/workspace/PPIDataProcess/GenerateOrtholog/src/org/blast/evalue/yeast_human/temp/human_geneid_new.txt'
    print 'clean and mask repeat id...'
    proteinMask(output, maskidout)
    print 'done'
    
    mapoutput='/home/yqian33/workspace/PPIDataProcess/GenerateOrtholog/src/org/blast/evalue/yeast_human/temp/human_gid_uac.txt'
    print 'map gene id of species2 to UniprotAC...'
    mapGid2Uac(maskidout,mapinput, mapoutput)
    print 'done'
    
    repeatoutput='/home/yqian33/workspace/PPIDataProcess/GenerateOrtholog/src/org/blast/evalue/yeast_human/temp/human_gid_uac_new.txt'
    print 'removing repeat map pairs...'
    repeatMask(mapoutput, repeatoutput)
    print 'done'
    
    ortout='/home/yqian33/workspace/PPIDataProcess/GenerateOrtholog/src/org/blast/evalue/yeast_human/yeast_human.ort'
    print 'get the final ortholog file identified by UniprotAC...'
    map(output1, ortout, repeatoutput)
    print 'Converstion finished!'
    

if __name__ == '__main__':
    idmapping()