'''
Created on Sep 3, 2014

@author: yqian33
'''
def i2d_dip():
    #read interactions in dip to a dict
    ins = open("/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/yeast.DIP.nif", 'r')
    pdict={}
    i=0
    for line in ins:
        items = line.split('\t')
        str=items[0]+items[1]
        pdict[str]=i
        i+=1
    ins.close()
    print i
    #read interaction from I2D and check one by one
    ous = open('/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/add/yeast.i2d.dip.nif', 'w')
    ous2 = open('/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/add/yeast.DIP.new.nif', 'w')
    
    ins2 = open("/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/yeast.i2d.nif", 'r')
    all=0
    count = 0
    count1= 0
    for line2 in ins2:
        all += 1
        items2 = line2.split('\t')
        str1=items2[0]+items2[1]
        str2=items2[1]+items2[0]
        if str1 in pdict or str2 in pdict:
            item = line2.split('\t')
            count1 += 1
            ous2.write(item[0]+'\t'+item[1]+'\t'+'1'+'\n')
            
        else:
            ous.write(line2)
            count = count + 1
    
    print all
    print count1
    print count
    ins2.close() 
    ous.close()
    ous2.close()

def i2d_dip_check():
    #read interactions in dip to a dict
    ins = open("/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/fly.DIP.nif", 'r')
    pdict={}
    i=0
    for line in ins:
        items = line.split('\t')
        str=items[0]+items[1]
        pdict[str]=i
        i+=1
    ins.close()
    print i
    #read interaction from I2D and check one by one
    #ous = open('/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/add/yeast.i2d.dip.nif', 'w')
    #ous2 = open('/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/add/yeast.DIP.new.nif', 'w')
    
    ins2 = open("/home/yqian33/APPI/software/AlignMCL-1.2/test/AlignMCLdataset/ppi/add/fly.DIP.new.nif", 'r')
    all=0
    count = 0
    count1= 0
    for line2 in ins2:
        all += 1
        items2 = line2.split('\t')
        str1=items2[0]+items2[1]
        str2=items2[1]+items2[0]
        if str1 in pdict or str2 in pdict:
            item = line2.split('\t')
            count1 += 1
            #ous2.write(item[0]+'\t'+item[1]+'\t'+'1'+'\n')
            
        else:
            #ous.write(line2)
            count = count + 1
    
    print all
    print count1
    print count
    ins2.close() 
    #ous.close()
    #ous2.close()
    
if __name__ == '__main__':
    i2d_dip_check()
    pass