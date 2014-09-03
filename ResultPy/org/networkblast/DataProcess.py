'''
Created on May 12, 2014

@author: yqian33
'''
def removeColumns():
    ins = open("/home/yqian33/APPI/ProSeqBla/idmapfile/YEAST_sysname_uniprotac.txt", 'r')
    ous = open("/home/yqian33/APPI/ProSeqBla/idmapfile/YEAST_sysname_uniprotac_NEW.txt", 'w')
    count = 0
    for line in ins:
        ous.write(line[75:])

    print count
    ins.close()
    ous.close()
    print 'done'
    
if __name__ == '__main__':
    removeColumns()
    pass