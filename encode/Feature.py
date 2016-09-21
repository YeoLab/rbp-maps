'''
Created on Sep 21, 2016

@author: brian
'''
import pybedtools as bt

class Feature():
    '''
    classdocs
    '''


    def __init__(self, line, source):
        '''
        Constructor
        '''
        self.line = line
        self.source = source

    def get_bedtool(self):
        if(self.source == 'miso'):
            chrom, start, end, strand = self.line.split(':')
            start = int(start) - 1
            end = int(end)
        return bt.create_interval_from_list([chrom,
                                             start,
                                             end,
                                             '0',
                                             '0',
                                             strand])

class SkippedExonFeature():
    def __init__(self, line, source):
        self.source = source
        self.line = line
    def get_bedtools(self):
        if(self.source == 'miso'):
            up, se, down = self.line.split('@')
            
            up = Feature(up,self.source).get_bedtool()
            se = Feature(se,self.source).get_bedtool()
            down = Feature(down,self.source).get_bedtool()
            
        return up, se, down

class A5ssFeature():
    def __init__(self, line, source):
        self.source = source
        self.line = line
    def get_bedtools(self):
        if(self.source == 'miso'):
            alt, downstream = self.line.split('@')
            chrom, start, end, strand = alt.split(':')
            end1, end2 = end.split('|')
            if(strand == '+'):
                splice1 = bt.create_interval_from_list([chrom,start,end1,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,end1,end2,'0','0',strand]) # middle
            else:
                splice1 = bt.create_interval_from_list([chrom,end2,start,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,end1,end2,'0','0',strand]) # middle
            chrom, start, end, strand = downstream.split(':')
            downstream = bt.create_interval_from_list([chrom,start,end,'0','0',strand])
        return splice1, splice2, downstream
    
class A3ssFeature():
    def __init__(self, line, source):
        self.source = source
        self.line = line
    def get_bedtools(self):
        if(self.source == 'miso'):
            upstream, alt = self.line.split('@')
            chrom, start, end, strand = upstream.split(':')
            upstream = bt.create_interval_from_list([chrom,start,end,'0','0',strand])
            
            chrom, start, end, strand = alt.split(':')
            start1, start2 = start.split('|')
            
            if(strand == '+'):
                splice1 = bt.create_interval_from_list([chrom,start1,start2,'0','0',strand]) # the middle one
                splice2 = bt.create_interval_from_list([chrom,start2,end,'0','0',strand]) # the downstream one
            elif(strand == '-'):
                splice1 = bt.create_interval_from_list([chrom,start1,start2,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,end,start1,'0','0',strand])
        return upstream, splice1, splice2
line = 'chr3:53274267:53274364:-@chr3:53271813:53271836:-@chr3:53268999:53269190:-'

F = SkippedExonFeature(line, 'miso')
up, se, down = F.get_bedtools()
print(up)
print(se)
print(down)

line = 'chr2:183800103:183799993|183800021:-@chr2:183799480:183799560:-'
F = A5ssFeature(line,'miso')
splice1, splice2, downstream = F.get_bedtools()
print(splice1)
print(splice2)
print(downstream)

line = 'chr2:55764619:55764721:+@chr2:55771074|55771161:55771210:+'
F = A3ssFeature(line,'miso')
upstream, splice1, splice2 = F.get_bedtools()
print(upstream)
print(splice1)
print(splice2)

