'''
Created on Jul 4, 2014

@author: yqian33
'''
#As the program for transfer PPI to Martin needed PPI can't have interaction scores
#so use this program to remove the score of each interaction
def nif_remove_score():
    ins = open('/home/yqian33/APPI/software/to_Yu/yiwei/neg.txt', 'r')
    ous = open('/home/yqian33/APPI/software/to_Yu/yiwei/worm.pairs.txt', 'w')
    
    for line in ins:
        item = line.split('\t')
        if item[1][-1] == '\n':
            ous.write(item[0]+'\t'+item[1])
        else:
            ous.write(item[0]+'\t'+item[1]+'\n')
        
    ins.close()
    ous.close()
    print 'done'
    
def parse_pred():
    ins = open('/home/yqian33/APPI/software/to_Yu/Martin_ppi_program/example/HELICO_preds.txt', 'r')
    ous = open('/home/yqian33/APPI/software/to_Yu/Martin_ppi_program/example/HELICO_preds_pairs.txt', 'w')
    
    linecount=0
    for line in ins:
        if linecount>1:
            item = line.split('\t')
            number = linecount-2
            for it in item[:-1]:
                str1 = str(linecount-2)
                str2 = str(number)
                ous.write(str1+"\t"+str2+"\t"+it+"\n")
                number = number + 1
        
        linecount = linecount + 1 
        
    ins.close()
    ous.close()
    print 'done'

def score_known_ppi():
    ins1 = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/yeast_predALL.txt', 'r')
    ins2 = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/yeast.DIP.nif', 'r')
    ous = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/yeast.DIP_score1.txt', 'w')
    
    dict = {}
    for line1 in ins1:
        items = line1.split("\t")
        pstr = items[0]+items[1]
        dict[pstr]=items[2][:-1]
    
    count=0
    for line in ins2:
        item = line.split("\t")
        #check item[0] and item[1]
        score='0'
        qstr=''
        if dict.has_key(item[0]+item[1]):
            score = dict[item[0]+item[1]]
            qstr=item[0]+'\t'+item[1]
        elif dict.has_key(item[1]+item[0]):
            score = dict[item[1]+item[0]]
            qstr=item[1]+'\t'+item[0]
        else:
            score = '0'
            qstr=item[0]+'\t'+item[1]
            count = count + 1
            
        ous.write(qstr+'\t'+score+'\n')
    
    ins1.close()
    ins2.close()
    ous.close()
    print count
    print 'done'

def testdataset_count():
    #ins = open('/home/yqian33/APPI/software/to_Yu/Martin_ppi_program/example/test.dat', 'r') 
    ins = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast.DIP/test.dat', 'r') 
    
    count=0
    count1=0
    for line in ins:
        if line[0] == '-':
            count = count + 1
        elif line[0] =='1':
            count1 = count1 + 1
    
    print '-:',count
    print '+:', count1
    ins.close()
    print 'done'
        
def generate_test_dataset():
    ins = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast.DIP/indicator.txt', 'r')   
    ous = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast.DIP/test.dat', 'w')   
    
    count=0
    for line in ins:
        if count==1:
            ous.write(line)
        count = count + 1
        if count == 9:
            count=0
            
    ins.close()   
    print 'done'   

class PPI:
    protein1=""
    protein2=""
    score=0.0
    def __init__(self, p1, p2, score):
        self.protein1=p1
        self.protein2=p2
        self.score = score
    def __del__(self):
        pass
    
    def getProtein1(self):
        return self.protein1
    def getProtein2(self):
        return self.protein2
    def getScore(self):
        return self.score

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
     
def compute_max_interaction():
    import operator
    
    infile = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/yeast.DIP.nif', 'r')
    dict={}
    i=0
    for pair in infile:
        item = pair.split('\t')
        p=item[0]+'\t'+item[1]
        dict[p]=i
        #print p
        i+=1
    infile.close()
        
    
    ins = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/yeast.DIP_score.txt', 'r')
    ous = open('/home/yqian33/workspace/PPIDataProcess/AddInteraction/src/yeast-test/yeast.DIP_score_sort.txt', 'w')
    
    count = 0
    ppiList=[]
    for line in ins:
        #print line
        item = line.split('\t')
        p=item[0]+'\t'+item[1]
        p1 = item[1]+'\t'+item[0]
        if is_number(item[2]):# and (p not in dict or p1 not in dict):
            ppi=PPI(item[0], item[1], float(item[2]))
            ppiList.append(ppi)
    
    ppiList.sort(key=operator.attrgetter('score'), reverse=True)
    
    for ppi1 in ppiList:
        ous.write(ppi1.getProtein1()+"\t"+ppi1.getProtein2()+"\t"+str(ppi1.getScore())+"\n")  
    
    ins.close()
    ous.close()
    print 'done'
    

if __name__ == '__main__':
    #parse_pred()
    #score_known_ppi()
    #nif_remove_score()
    #testdataset_count()
    #generate_test_dataset()
    #testdataset_count()
    compute_max_interaction()
    pass