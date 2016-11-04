'''
Created on Sep 21, 2016

@author: brian
'''
import pybedtools as bt

class Feature():
    '''
    classdocs
    '''

    def __init__(self, annotation, source):
        '''
        Constructor
        '''
        self.source = source
        self.annotation = annotation.rstrip()
        
    def get_bedtool(self):
        if(self.source == 'bed'):
            chrom, start, end, name, score, strand = self.annotation.split('\t')
        return bt.create_interval_from_list([chrom,
                                             start,
                                             end,
                                             name,
                                             score,
                                             strand])

class SkippedExonFeature(Feature):
    def __init__(self, annotation, source):
        Feature.__init__(self, annotation, source)

    def get_bedtools(self):
        if(self.source == 'miso'):
            event = self.annotation.split('\t')[0]
            up, se, down = event.split('@')
            
            chrom, start, stop, strand = up.split(':')
            up = bt.create_interval_from_list([chrom, int(start)-1, stop, '0', '0', strand])
            
            chrom, start, stop, strand = se.split(':')
            se = bt.create_interval_from_list([chrom, int(start)-1, stop, '0', '0', strand])
            
            chrom, start, stop, strand = down.split(':')
            down = bt.create_interval_from_list([chrom, int(start)-1, stop, '0', '0', strand])
        elif(self.source == 'hta2_0'):
            pass
        elif(self.source == 'xintao'):
            pass
        elif(self.source == 'eric'):
            name, se = self.annotation.split(';')
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
        elif(self.source == 'rmats'):
            id, GeneID, geneSymbol, chrom, strand, \
            exonStart_0base, exonEnd, \
            upstreamES, upstreamEE, \
            downstreamES, downstreamEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, \
            IncFormLen, SkipFormLen, PValue, \
            FDR, IncLevel1, IncLevel2, IncLevelDifference = self.annotation.split('\t')
            
            se = bt.create_interval_from_list([chrom, exonStart_0base, exonEnd, '0', '0', strand])
            if(strand == '+'):
                up = bt.create_interval_from_list([chrom, upstreamES, upstreamEE, '0', '0', strand])
                down = bt.create_interval_from_list([chrom, downstreamES, downstreamEE, '0', '0', strand])
            elif(strand == '-'):
                down = bt.create_interval_from_list([chrom, upstreamES, upstreamEE, '0', '0', strand])
                up = bt.create_interval_from_list([chrom, downstreamES, downstreamEE, '0', '0', strand])
            else:
                print("Warning, strand not correct!")
                return -1
        return up, se, down

class A5ssFeature():
    def __init__(self, annotation, source):
        Feature.__init__(self, annotation, source)

    def get_bedtools(self):
        if(self.source == 'miso'):
            event = self.annotation.split('\t')[0]
            alt, downstream = event.split('@')
            chrom, start, end, strand = alt.split(':')
            end1, end2 = end.split('|')
            if(strand == '+'):
                splice1 = bt.create_interval_from_list([chrom,int(start)-1,end1,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,int(start)-1,end2,'0','0',strand]) # middle
            else:
                splice1 = bt.create_interval_from_list([chrom,int(end2)-1,start,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,int(end1)-1,start,'0','0',strand]) # middle
            chrom, start, end, strand = downstream.split(':')
            downstream = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
        elif(self.source == 'rmats'):
            ID, GeneID, geneSymbol, chrom, strand, \
            longExonStart_0base, longExonEnd, \
            shortES, shortEE, \
            flankingES, flankingEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, \
            IncFormLen, SkipFormLen, \
            PValue, FDR, \
            IncLevel1, IncLevel2, IncLevelDifference = self.annotation.split('\t')
            
            downstream = bt.create_interval_from_list([chrom, flankingES, flankingEE, '0', '0', strand])
            splice1 = bt.create_interval_from_list([chrom, longExonStart_0base, longExonEnd, '0', '0', strand])
            splice2 = bt.create_interval_from_list([chrom, shortES, shortEE, '0', '0', strand])
            
        return splice1, splice2, downstream
    
class A3ssFeature():
    def __init__(self, annotation, source):
        Feature.__init__(self, annotation, source)

    def get_bedtools(self):
        """
        Produces 3 intervals: 
        - upstream exon
        - splice1 (the 'longer' one)
        - splice2 (the 'downstream' one)
        
        eg: 
        
        [upstream exon]------------[         splice1] (+)
        [upstream exon]---------------------[splice2] (+)
        
        [splice1         ]------------[upstream exon] (-)
        [splice2]---------------------[upstream exon] (-)
        """
        if(self.source == 'miso'):
            event = self.annotation.split('\t')[0]
            upstream, alt = event.split('@')
            chrom, start, end, strand = upstream.split(':')
            upstream = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
            
            chrom, start, end, strand = alt.split(':')
            start1, start2 = start.split('|')
            
            if(strand == '+'):
                splice1 = bt.create_interval_from_list([chrom,int(start1)-1,end,'0','0',strand]) # the longer one
                splice2 = bt.create_interval_from_list([chrom,int(start2)-1,end,'0','0',strand]) # the shorter one
            elif(strand == '-'):
                splice1 = bt.create_interval_from_list([chrom,int(end)-1,start2,'0','0',strand])
                splice2 = bt.create_interval_from_list([chrom,int(end)-1,start1,'0','0',strand])
        elif(self.source == 'rmats'):
            ID, GeneID, geneSymbol, chrom, strand, \
            longExonStart_0base, longExonEnd, \
            shortES, shortEE, \
            flankingES, flankingEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, IncFormLen, SkipFormLen, \
            PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference = self.annotation.split('\t')
            
            upstream = bt.create_interval_from_list([chrom, flankingES, flankingEE, '0', '0', strand])
            splice1 = bt.create_interval_from_list([chrom, longExonStart_0base, longExonEnd, '0', '0', strand])
            splice2 = bt.create_interval_from_list([chrom, shortES, shortEE, '0', '0', strand])
        return upstream, splice1, splice2


class RIFeature():
    def __init__(self, annotation, source):
        Feature.__init__(self, annotation, source)

    def get_bedtools(self):
        
        if(self.source == 'xintao'):
            """
            I THINK THESE ARE ZERO-BASED, BUT I'M NOT SURE...
            
            CCT8_ENSG00000156261.8;RI:chr21:30434649:30434736-30434811:30434896:-
            EXOSC8_ENSG00000120699.8;RI:chr13:37577071:37577144-37578614:37578698:+
            """
            annotation, chrom, region1, region2, region3, strand = self.annotation.rstrip().split(':')
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
        elif(self.source == 'rmats'):
            ID, GeneID, geneSymbol, chrom, strand, \
            riExonStart_0base, riExonEnd, \
            upstreamES, upstreamEE, \
            downstreamES, downstreamEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, \
            IncFormLen, SkipFormLen, \
            PValue, FDR, \
            IncLevel1, IncLevel2, IncLevelDifference = self.annotation.split('\t')
            if(strand == '+'):
                upstream = bt.create_interval_from_list([chrom, upstreamES, upstreamEE, '0', '0', strand])
                downstream = bt.create_interval_from_list([chrom, downstreamES, downstreamEE, '0', '0', strand])
            elif(strand == '-'):
                downstream = bt.create_interval_from_list([chrom, upstreamES, upstreamEE, '0', '0', strand])
                upstream = bt.create_interval_from_list([chrom, downstreamES, downstreamEE, '0', '0', strand])
            else:
                print("strand not correct")
                return -1
        return upstream, downstream
class MXEFeature():
    def __init__(self, annotation, source):
        Feature.__init__(self, annotation, source)

    def get_bedtools(self):
        
        if(self.source == 'rmats'):
            """
            1stExon describes the upstream mutually exclusive exon
            2ndExon describes the downstream mutually excluxive exon
            """
            ID, GeneID, geneSymbol, chrom, strand, \
            firstExonStart_0base, firstExonEnd, \
            secondExonStart_0base, secondExonEnd, \
            upstreamES, upstreamEE, \
            downstreamES, downstreamEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, IncFormLen, SkipFormLen, \
            PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference = self.annotation.split('\t')
            
            if(strand == '+'):
                upstream = bt.create_interval_from_list([chrom, upstreamES, upstreamEE, '0', '0', strand])
                downstream = bt.create_interval_from_list([chrom, downstreamES, downstreamEE, '0', '0', strand])
                up_mxe = bt.create_interval_from_list([chrom, firstExonStart_0base, firstExonEnd, '0', '0', strand])
                down_mxe = bt.create_interval_from_list([chrom, secondExonStart_0base, secondExonEnd, '0', '0', strand])
            elif(strand == '-'): # upstream/downstream is flipped for rmats
                downstream = bt.create_interval_from_list([chrom, upstreamES, upstreamEE, '0', '0', strand])
                upstream = bt.create_interval_from_list([chrom, downstreamES, downstreamEE, '0', '0', strand])
                down_mxe = bt.create_interval_from_list([chrom, firstExonStart_0base, firstExonEnd, '0', '0', strand])
                up_mxe = bt.create_interval_from_list([chrom, secondExonStart_0base, secondExonEnd, '0', '0', strand])
            else:
                print("Warning, strand not correct!")
                return -1
        return upstream, up_mxe, down_mxe, downstream
"""        
annotation = 'chr3:53274267:53274364:-@chr3:53271813:53271836:-@chr3:53268999:53269190:-'
print(annotation)
F = SkippedExonFeature(annotation, 'miso')
up, se, down = F.get_bedtools()
print(up)
print(se)
print(down)

annotation = 'chr10:102743831:102743705|102743791:-@chr10:102743512:102743574:-'
print(annotation)
F = A5ssFeature(annotation,'miso')
splice1, splice2, downstream = F.get_bedtools()
print(splice1)
print(splice2)
print(downstream)

annotation = 'chr17:80417868:80417948|80418199:+@chr17:80422163:80422306:+'
print(annotation)
F = A5ssFeature(annotation,'miso')
splice1, splice2, downstream = F.get_bedtools()
print(splice1)
print(splice2)
print(downstream)

annotation = 'CCT8_ENSG00000156261.8;RI:chr21:30434649:30434736-30434811:30434896:-'
print(annotation)
F = RIFeature(annotation,'xintao')
upstream, downstream = F.get_bedtools()
print(upstream)
print(downstream)
annotation = 'EXOSC8_ENSG00000120699.8;RI:chr13:37577071:37577144-37578614:37578698:+'
print(annotation)
F = RIFeature(annotation,'xintao')
upstream, downstream = F.get_bedtools()
print(upstream)
print(downstream)

annotation = 'chr10:100185575:100185742:-@chr10:100185441|100185477:100185298:-'
print(annotation)
F = A3ssFeature(annotation,'miso')
upstream, splice1, splice2 = F.get_bedtools()
print(upstream)
print(splice1)
print(splice2)

annotation = 'chr2:55764619:55764721:+@chr2:55771074|55771161:55771210:+'
print(annotation)
F = A3ssFeature(annotation,'miso')
upstream, splice1, splice2 = F.get_bedtools()
print(upstream)
print(splice1)
print(splice2)

annotation = 'GPBP1_ENSG00000062194.11;A3SS:chr5:56543042-56545236:56543042-56545281:+||ENST00000513524.1:chr5:56542902-56543042:+||ENST00000511209.1:chr5:56545281-56545403:+'
print(annotation)
F = A3ssFeature(annotation,'eric')
upstream, splice1, splice2 = F.get_bedtools()
print(upstream)
print(splice1)
print(splice2)

annotation = 'SLC26A6_ENSG00000225697.6;A3SS:chr3:48665379-48664540:48665379:48664485:-||Not_found||Not_found'
print(annotation)
F = A3ssFeature(annotation,'eric')
upstream, splice1, splice2 = F.get_bedtools()
print(upstream)
print(splice1)
print(splice2)
"""