#!/bin/env python

'''
Created on Sep 21, 2016

@author: brian
'''
import pybedtools as bt


class Feature():
    '''
    classdocs
    '''

    def __init__(self, annotation_line, annotation_format):
        '''
        Constructor
        '''
        self.source = annotation_format
        self.annotation = annotation_line.rstrip()

    def get_bedtool(self):
        chrom = None
        start = None
        end = None
        name = None
        score = None
        strand = None
        if self.source == 'bed':
            chrom, start, end, name, score, strand = \
                self.annotation.split('\t')
        return bt.create_interval_from_list(
            [chrom, start, end, name, score, strand]
        )


class Skipped_exon(Feature):

    def __init__(self, annotation_line, annotation_format):
        Feature.__init__(self, annotation_line, annotation_format)

    def get_bedtools(self):
        """
        Returns an upstream, skipped, and downstream exon interval.

        Returns
        -------

        """
        up = None
        se = None
        down = None

        if self.source == 'miso':
            event = self.annotation.split('\t')[0]
            up, se, down = event.split('@')

            chrom, start, stop, strand = up.split(':')
            up = bt.create_interval_from_list(
                [chrom, int(start) - 1, stop, '0', '0', strand]
            )

            chrom, start, stop, strand = se.split(':')
            se = bt.create_interval_from_list(
                [chrom, int(start) - 1, stop, '0', '0', strand]
            )

            chrom, start, stop, strand = down.split(':')
            down = bt.create_interval_from_list(
                [chrom, int(start) - 1, stop, '0', '0', strand]
            )
        elif self.source == 'hta2_0':
            pass
        elif self.source == 'xintao':
            pass
        elif self.source == 'eric':
            pos, upstream, skipped, downstream = self.annotation.split('\t')
            chrom, strand, up_junc, down_junc, se_region = pos.split('|')

            if strand == '-':  # these are coord-based not *stream-based.
                down_start, down_end = upstream.split('-')
                up_start, up_end = downstream.split('-')
            else:
                up_start, up_end = upstream.split('-')
                down_start, down_end = downstream.split('-')
            se_start, se_end = skipped.split('-')

            up = bt.create_interval_from_list(
                [chrom, up_start, up_end, '0', '0', strand]
            )
            down = bt.create_interval_from_list(
                [chrom, down_start, down_end, '0', '0', strand]
            )
            se = bt.create_interval_from_list(
                [chrom, se_start, se_end, '0', '0', strand]
            )
        elif self.source == 'rmats':
            # TODO: replace unused variables with _
            id, GeneID, geneSymbol, chrom, strand, \
            exonStart_0base, exonEnd, \
            upstreamES, upstreamEE, \
            downstreamES, downstreamEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, \
            IncFormLen, SkipFormLen, PValue, \
            FDR, IncLevel1, IncLevel2, IncLevelDifference = \
                self.annotation.split('\t')

            se = bt.create_interval_from_list(
                [chrom, exonStart_0base, exonEnd, '0', '0', strand]
            )
            if strand == '+':
                up = bt.create_interval_from_list(
                    [chrom, upstreamES, upstreamEE, '0', '0', strand]
                )
                down = bt.create_interval_from_list(
                    [chrom, downstreamES, downstreamEE, '0', '0', strand]
                )
            elif strand == '-':
                down = bt.create_interval_from_list(
                    [chrom, upstreamES, upstreamEE, '0', '0', strand]
                )
                up = bt.create_interval_from_list(
                    [chrom, downstreamES, downstreamEE, '0', '0', strand]
                )
            else:
                print("Warning, strand not correct!")
                return 1
        return up, se, down


class Alt_5p_splice_site(Feature):
    def __init__(self, annotation_line, annotation_format):
        Feature.__init__(self, annotation_line, annotation_format)

    def get_bedtools(self):
        """
        Produces 3 intervals:
        - splice1 (the 'longer' one)
        - splice2 (the 'downstream' one)
        - downstream exon

        eg:

        [splice1         ]------------[downstream exon] (+)
        [splice2]---------------------[downstream exon] (+)

        [downstream exon]------------[         splice1] (-)
        [downstream exon]---------------------[splice2] (-)

        Returns
        -------

        splice1 : 
        splice2 :
        downstream :
        """
        splice1 = None
        splice2 = None
        downstream = None

        if self.source == 'miso':
            event = self.annotation.split('\t')[0]
            alt, downstream = event.split('@')
            chrom, start, end, strand = alt.split(':')
            end1, end2 = end.split('|')
            if strand == '+':
                splice1 = bt.create_interval_from_list(
                    [chrom, int(start) - 1, end1, '0', '0', strand]
                )
                splice2 = bt.create_interval_from_list(
                    [chrom, int(start) - 1, end2, '0', '0', strand]
                )  # middle
            else:
                splice1 = bt.create_interval_from_list(
                    [chrom, int(end2) - 1, start, '0', '0', strand]
                )
                splice2 = bt.create_interval_from_list(
                    [chrom, int(end1) - 1, start, '0', '0', strand]
                )  # middle
            chrom, start, end, strand = downstream.split(':')
            downstream = bt.create_interval_from_list(
                [chrom, int(start) - 1, end, '0', '0', strand]
            )
        elif self.source == 'rmats':
            ID, GeneID, geneSymbol, chrom, strand, \
            longExonStart_0base, longExonEnd, \
            shortES, shortEE, \
            flankingES, flankingEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, \
            IncFormLen, SkipFormLen, \
            PValue, FDR, \
            IncLevel1, IncLevel2, IncLevelDifference = \
                self.annotation.split('\t')

            downstream = bt.create_interval_from_list(
                [chrom, flankingES, flankingEE, '0', '0', strand]
            )
            splice1 = bt.create_interval_from_list(
                [chrom, longExonStart_0base, longExonEnd, '0', '0', strand]
            )
            splice2 = bt.create_interval_from_list(
                [chrom, shortES, shortEE, '0', '0', strand]
            )

        return splice1, splice2, downstream


class Alt_3p_splice_site(Feature):
    def __init__(self, annotation_line, annotation_format):
        Feature.__init__(self, annotation_line, annotation_format)

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
        upstream = None
        splice1 = None
        splice2 = None

        if self.source == 'miso':
            event = self.annotation.split('\t')[0]
            upstream, alt = event.split('@')
            chrom, start, end, strand = upstream.split(':')
            upstream = bt.create_interval_from_list(
                [chrom, int(start) - 1, end, '0', '0', strand]
            )

            chrom, start, end, strand = alt.split(':')
            start1, start2 = start.split('|')

            if strand == '+':
                splice1 = bt.create_interval_from_list(
                    [chrom, int(start1) - 1, end, '0', '0', strand]
                )  # the longer one
                splice2 = bt.create_interval_from_list(
                    [chrom, int(start2) - 1, end, '0', '0', strand]
                )  # the shorter one
            elif strand == '-':
                splice1 = bt.create_interval_from_list(
                    [chrom, int(end) - 1, start2, '0', '0', strand]
                )
                splice2 = bt.create_interval_from_list(
                    [chrom, int(end) - 1, start1, '0', '0', strand]
                )
        elif self.source == 'rmats':
            ID, GeneID, geneSymbol, chrom, strand, \
            longExonStart_0base, longExonEnd, \
            shortES, shortEE, \
            flankingES, flankingEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, IncFormLen, SkipFormLen, \
            PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference = \
                self.annotation.split('\t')

            upstream = bt.create_interval_from_list(
                [chrom, flankingES, flankingEE, '0', '0', strand]
            )
            splice1 = bt.create_interval_from_list(
                [chrom, longExonStart_0base, longExonEnd, '0', '0', strand]
            )
            splice2 = bt.create_interval_from_list(
                [chrom, shortES, shortEE, '0', '0', strand]
            )
        return upstream, splice1, splice2


class Retained_intron(Feature):
    def __init__(self, annotation_line, annotation_format):
        Feature.__init__(self, annotation_line, annotation_format)

    def get_bedtools(self):
        upstream = None
        downstream = None

        if self.source == 'xintao':
            """
            I THINK THESE ARE ZERO-BASED, BUT I'M NOT SURE...
            
            CCT8_ENSG00000156261.8;RI:chr21:30434649:30434736-30434811:30434896:-
            EXOSC8_ENSG00000120699.8;RI:chr13:37577071:37577144-37578614:37578698:+
            """
            annotation, chrom, region1, region2, region3, strand = \
                self.annotation.rstrip().split(':')
            if strand == '+':
                upstream_start = region1
                upstream_end, downstream_start = region2.split('-')
                downstream_end = region3
            elif strand == '-':
                downstream_start = region1
                downstream_end, upstream_start = region2.split('-')
                upstream_end = region3
            else:
                print("invalid strand information, defaulting to +")
                upstream_start = region1
                upstream_end, downstream_start = region2.split('-')
                downstream_end = region3
            upstream = bt.create_interval_from_list(
                [chrom, upstream_start, upstream_end, '0', '0', strand]
            )
            downstream = bt.create_interval_from_list(
                [chrom, downstream_start, downstream_end, '0', '0', strand]
            )
        elif self.source == 'rmats':
            ID, GeneID, geneSymbol, chrom, strand, \
            riExonStart_0base, riExonEnd, \
            upstreamES, upstreamEE, \
            downstreamES, downstreamEE, \
            ID1, IJC_SAMPLE_1, SJC_SAMPLE_1, \
            IJC_SAMPLE_2, SJC_SAMPLE_2, \
            IncFormLen, SkipFormLen, \
            PValue, FDR, \
            IncLevel1, IncLevel2, IncLevelDifference = \
                self.annotation.split('\t')
            if strand == '+':
                upstream = bt.create_interval_from_list(
                    [chrom, upstreamES, upstreamEE, '0', '0', strand]
                )
                downstream = bt.create_interval_from_list(
                    [chrom, downstreamES, downstreamEE, '0', '0', strand]
                )
            elif strand == '-':
                downstream = bt.create_interval_from_list(
                    [chrom, upstreamES, upstreamEE, '0', '0', strand]
                )
                upstream = bt.create_interval_from_list(
                    [chrom, downstreamES, downstreamEE, '0', '0', strand]
                )
            else:
                print("strand not correct")
                return -1
        return upstream, downstream


class Mutually_exclusive_exon(Feature):
    def __init__(self, annotation_line, annotation_format):
        Feature.__init__(self, annotation_line, annotation_format)

    def get_bedtools(self):
        upstream = None
        up_mxe = None
        down_mxe = None
        downstream = None

        if self.source == 'rmats':
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
            PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference = \
                self.annotation.split('\t')

            if strand == '+':
                upstream = bt.create_interval_from_list(
                    [chrom, upstreamES, upstreamEE, '0', '0', strand])
                downstream = bt.create_interval_from_list(
                    [chrom, downstreamES, downstreamEE, '0', '0', strand])
                up_mxe = bt.create_interval_from_list(
                    [chrom, firstExonStart_0base, firstExonEnd, '0', '0',
                     strand]
                )
                down_mxe = bt.create_interval_from_list(
                    [chrom, secondExonStart_0base, secondExonEnd, '0', '0',
                     strand]
                )
            elif strand == '-':  # upstream/downstream is flipped for rmats
                downstream = bt.create_interval_from_list(
                    [chrom, upstreamES, upstreamEE, '0', '0', strand])
                upstream = bt.create_interval_from_list(
                    [chrom, downstreamES, downstreamEE, '0', '0', strand])
                down_mxe = bt.create_interval_from_list(
                    [chrom, firstExonStart_0base, firstExonEnd, '0', '0',
                     strand]
                )
                up_mxe = bt.create_interval_from_list(
                    [chrom, secondExonStart_0base, secondExonEnd, '0', '0',
                     strand]
                )
            else:
                print("Warning, strand not correct!")
                return -1
        return upstream, up_mxe, down_mxe, downstream


class UnscaledCDS(Feature):
    def __init__(self, annotation_line, annotation_format):
        Feature.__init__(self, annotation_line, annotation_format)

    def get_bedtools(self):
        """
        This is a 12-column, tabbed line basically concatenating two BED6
        features.
        Parameters
        ----------

        Returns
        -------
        upstream : pybedtools.BedTool
        downstream : pybedtools.BedTool

        """
        if self.source == 'twobeds':
            first_chrom, first_start, first_end, \
            first_name, first_score, first_strand, \
            second_chrom, second_start, second_end, \
            second_name, second_score, second_strand = \
                self.annotation.split('\t')
            if first_strand == '+':
                upstream = bt.create_interval_from_list(
                    [first_chrom, first_start, first_end,
                     first_name, first_score, first_strand]
                )
                downstream = bt.create_interval_from_list(
                    [second_chrom, second_start, second_end,
                     second_name, second_score, second_strand]
                )
            else:
                downstream = bt.create_interval_from_list(
                    [first_chrom, first_start, first_end,
                     first_name, first_score, first_strand]
                )
                upstream = bt.create_interval_from_list(
                    [second_chrom, second_start, second_end,
                     second_name, second_score, second_strand]
                )
        else:
            print("Warning, no valid feature source found. ")
        return upstream, downstream