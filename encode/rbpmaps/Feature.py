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
        self.line = line.rstrip()
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
        self.line = line.rstrip()
    def get_bedtools(self):
        if(self.source == 'miso'):
            up, se, down = self.line.split('@')
            
            chrom, start, stop, strand = up.split(':')
            up = bt.create_interval_from_list([chrom, int(start)-1, stop, '0', '0', strand])
            
            chrom, start, stop, strand = se.split(':')
            se = bt.create_interval_from_list([chrom, int(start)-1, stop, '0', '0', strand])
            
            chrom, start, stop, strand = down.split(':')
            down = bt.create_interval_from_list([chrom, int(start)-1, stop, '0', '0', strand])
        elif(self.source == 'xintao'):
            pass
        elif(self.source == 'eric'):
            name, se = self.line.split(';')
            xintao, ericleft, ericright = se.split('||')
            upstream_es = 1
            downstream_es = 250000000
            if("Not_found") not in ericleft:
                upstream_es = ericleft.split(':')[2].split('-')[0]
            if("Not_found") not in ericright:
                downstream_ee = ericright.split(':')[2].split('-')[1]
            
            event, chrom, upstream, downstream, strand = xintao.split(':')
            upstream_ee, skipped_es = upstream.split('-')
            skipped_ee, downstream_es = downstream.split('-')
            
            se = bt.create_interval_from_list([chrom, skipped_es, skipped_ee, '0', '0', strand])
            if(strand == '+'):
                up = bt.create_interval_from_list([chrom, upstream_es, upstream_ee, '0', '0', strand])
                down = bt.create_interval_from_list([chrom, downstream_es, downstream_ee, '0', '0', strand])
            elif(strand == '-'):
                up = bt.create_interval_from_list([chrom, downstream_es, downstream_ee, '0', '0', strand])
                down = bt.create_interval_from_list([chrom, upstream_es, upstream_ee, '0', '0', strand])
        return up, se, down

class A5ssFeature():
    def __init__(self, line, source):
        self.source = source
        self.line = line.rstrip()
    def get_bedtools(self):
        if(self.source == 'miso'):
            alt, downstream = self.line.split('@')
            chrom, start, end, strand = alt.split(':')
            end1, end2 = end.split('|')
            if(strand == '+'):
                splice1 = bt.create_interval_from_list([chrom,int(start)-1,end1,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,int(end1)-1,end2,'0','0',strand]) # middle
            else:
                splice1 = bt.create_interval_from_list([chrom,int(end2)-1,start,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,int(end1)-1,end2,'0','0',strand]) # middle
            chrom, start, end, strand = downstream.split(':')
            downstream = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
        return splice1, splice2, downstream
    
class A3ssFeature():
    def __init__(self, line, source):
        self.source = source
        self.line = line.rstrip()
    def get_bedtools(self):
        if(self.source == 'miso'):
            upstream, alt = self.line.split('@')
            chrom, start, end, strand = upstream.split(':')
            upstream = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
            
            chrom, start, end, strand = alt.split(':')
            start1, start2 = start.split('|')
            
            if(strand == '+'):
                splice1 = bt.create_interval_from_list([chrom,int(start1)-1,start2,'0','0',strand]) # the middle one
                splice2 = bt.create_interval_from_list([chrom,int(start2)-1,end,'0','0',strand]) # the downstream one
            elif(strand == '-'):
                splice1 = bt.create_interval_from_list([chrom,int(start1)-1,start2,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,int(end)-1,start1,'0','0',strand])
        return upstream, splice1, splice2

class RIFeature():
    def __init__(self, line, source):
        self.source = source
        self.line = line.rstrip()
    def get_bedtools(self):
        
        if(self.source == 'xintao'):
            """
            CCT8_ENSG00000156261.8;RI:chr21:30434649:30434736-30434811:30434896:-
            EXOSC8_ENSG00000120699.8;RI:chr13:37577071:37577144-37578614:37578698:+
            """
            annotation, chrom, region1, region2, region3, strand = self.line.rstrip().split(':')
            if(strand == '+'):
                upstream_start = region1
                upstream_end, downstream_start = region2.split('-')
                downstream_end = region3
            elif(strand == '-'):
                downstream_start = region1
                downstream_end, upstream_start = region2.split('-')
                upstream_end = region3
            else:
                print("invalid strand information, defaulting to +")
                upstream_start = region1
                upstream_end, downstream_start = region2.split('-')
                downstream_end = region3
            upstream = bt.create_interval_from_list([chrom,upstream_start,upstream_end,'0','0',strand])
            downstream = bt.create_interval_from_list([chrom,downstream_start,downstream_end,'0','0',strand]) 
        return upstream, downstream
        
"""line = 'chr3:53274267:53274364:-@chr3:53271813:53271836:-@chr3:53268999:53269190:-'

F = SkippedExonFeature(line, 'miso')
up, se, down = F.get_bedtools()
print(up)
print(se)
print(down)

line = 'chr10:102743831:102743705|102743791:-@chr10:102743512:102743574:-'
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

line = 'CCT8_ENSG00000156261.8;RI:chr21:30434649:30434736-30434811:30434896:-'
F = RIFeature(line,'xintao')
upstream, downstream = F.get_bedtools()
print(upstream)
print(downstream)"""

