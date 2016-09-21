'''
Created on Jun 18, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import pybedtools
import intervals
import misc
import Feature

def create_matrix(annotation, density, left = 300, right = 300, is_scaled = True):
    print("creating the matrix for {}".format(density.get_name()))
    print(left)
    print(right)
    # print("is this going to be scaled? {}".format(is_scaled))
    count = 0
    densities = {}
    # if(type(annotation) == pd.DataFrame)
    if(type(annotation) != pybedtools.bedtool.BedTool):
        bed_tool = misc.create_bedtool(annotation)
    else:
        bed_tool = annotation

    for interval in bed_tool:
        try:
            # print(interval)
            count = count + 1
            if count % 50000 == 0:
                print('processed {} features'.format(count))
            wiggle = pd.Series(intervals.some_range(density, interval, 0, 0))
            """
            We should process the intervals and the flanking regions separately 
            if we only want to scale the feature of interest. 
            """
            if not all(np.isnan(wiggle)):
                
                wiggle = wiggle.fillna(0) # convert all nans to 0
                
                wiggle = abs(wiggle) # convert all values to positive
                
                if(is_scaled == True):
                    wiggle = intervals.get_scale(wiggle)
                
                """
                For positive stranded features:
                    upstream interval is effectively 'left' of the interval,
                    downstream interval is effectively 'right', coordinate-wise.
                For negative stranded features:
                    upstream interval is effectively 'right' (coordinates are larger),
                    downstream interval is 'left' 
                """
                if interval.strand == '+': 
                    upstream_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                              interval.start-left,
                                                                              interval.start,
                                                                              interval.name,
                                                                              interval.score,
                                                                              interval.strand])
                    downstream_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                                interval.stop,
                                                                                interval.stop+right,
                                                                                interval.name,
                                                                                interval.score,
                                                                                interval.strand])
                else:
                    upstream_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                              interval.stop,
                                                                              interval.stop+right,
                                                                              interval.name,
                                                                              interval.score,
                                                                              interval.strand])
                    downstream_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                                interval.start-left,
                                                                                interval.start,
                                                                                interval.name,
                                                                                interval.score,
                                                                                interval.strand])
                
                upstream_wiggle = pd.Series(intervals.some_range(density,upstream_interval, 0, 0))
                upstream_wiggle = pd.Series(abs(np.nan_to_num(upstream_wiggle)))
                
                # print(type(upstream_wiggle))
                # print(type(wiggle))
                downstream_wiggle = pd.Series(intervals.some_range(density,downstream_interval, 0, 0))
                downstream_wiggle = pd.Series(abs(np.nan_to_num(downstream_wiggle)))
                
                # print(pd.concat([upstream_wiggle,wiggle,downstream_wiggle]).reset_index(drop=True).shape)
                # print(len(pd.concat([upstream_wiggle,wiggle,downstream_wiggle],axis=1).reset_index(drop=True)))
                densities[intervals.rename_index(interval)] = pd.concat([upstream_wiggle,wiggle,downstream_wiggle]).reset_index(drop=True)
                
                # densities[intervals.rename_index(interval)] = wiggle

        except Exception as e:
            print(e)
            print(type(wiggle))
            print(type(upstream_wiggle))
            print(type(downstream_wiggle))
    print("SUCCESS")
    return pd.DataFrame(densities).T

def create_a5ss_matrix(annotation, density, exon_offset, intron_offset, is_scaled, combine_regions = True):
    # chr17:80009218:80008888|80009170:-@chr17:80008538:80008640:-    ENSG00000169733
    # chr17:80417868:80417948|80418199:+@chr17:80422163:80422306:+    ENSG00000141562
    # chr2:55764619:55764721:+@chr2:55771074|55771161:55771210:+      ENSG00000163001
    # chr17:62502194:62502407:-@chr17:62500960|62500998:62500795:-    ENSG00000108654
    """
    Four regions:
    """
    # [    |    ]----|---|----[    |    [    |    ]
    three_alt1 = {}
    five_alt2 = {}
    three_alt2 = {}
    five_downstream = {}
    # Feature format:
    # chr17:80009218:80008888|80009170:-@chr17:80008538:80008640:-    ENSG00000169733
    # chr2:183800103:183799993|183800021:-@chr2:183799480:183799560:-    ENSG00000061676
    # chr17:80417868:80417948|80418199:+@chr17:80422163:80422306:+    ENSG00000141562
    with open(annotation) as f:
        # f.next() # for title
        for line in f:
            if not line.startswith('event_name'):
                event = line.split('\t')[0]
                alt1, alt2, downstream = Feature.A5ssFeature(event,'miso').get_bedtools()
                """three prime alt1 region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        alt1,
                                                                        alt2,
                                                                        exon_offset,
                                                                        0)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
        
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) 
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_alt1[event] = wiggle
                    
                """five prime of alt2 region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        alt1,
                                                                        alt2,
                                                                        exon_offset,
                                                                        0)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle)
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_alt2[event] = wiggle
                """three prime site of alt2 region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                         downstream,
                                                                         alt2,
                                                                         exon_offset,
                                                                         intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) #
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_alt2[event] = wiggle
                    # print("length of 3p skipped: {}".format(len(wiggle)))
        
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        alt2,
                                                                        downstream,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_downstream[event] = wiggle
                    # print("length of downstream: {}".format(len(wiggle)))
        three_alt1 = pd.DataFrame(three_alt1).T
        five_alt2 = pd.DataFrame(five_alt2).T
        three_alt2 = pd.DataFrame(three_alt2).T
        five_downstream = pd.DataFrame(five_downstream).T
        
        if combine_regions == False:
            return three_alt1, five_alt2, three_alt2, five_downstream
        else:
            ra = pd.concat([three_alt1,five_alt2,three_alt2,five_downstream],axis=1)
            ra.columns = range(0,ra.shape[1])
            return ra      
def create_a3ss_matrix(annotation, density, exon_offset, intron_offset, is_scaled, combine_regions=True):
    # chr2:55764619:55764721:+@chr2:55771074|55771161:55771210:+      ENSG00000163001
    # chr17:62502194:62502407:-@chr17:62500960|62500998:62500795:-    ENSG00000108654
    """
    Four regions:
    """
    # [    |    ]----|---|----[    |    [    |    ]
    three_upstream = {}
    five_alt1 = {}
    three_alt1 = {}
    five_alt2 = {}
    
    
    with open(annotation) as f:
        # f.next() # for title
        for line in f:
            if not line.startswith('event_name'):
                event = line.split('\t')[0]
                """ Deprecated 
                
                region1, region2 = event.split('@')
                
                upstream = misc.create_bed_tool_from_miso_a3ss(region1, False)
                alt1, alt2 = misc.create_bed_tool_from_miso_a3ss(region2, True)
                """
                upstream, alt1, alt2 = Feature.A3ssFeature(event,'miso').get_bedtools()
                # print('three prime site upstream')
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        alt1,
                                                                        upstream,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
        
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) 
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_upstream[event] = wiggle
        
                """five prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream,
                                                                        alt1,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle)
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_alt1[event] = wiggle
                
                """three prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                         alt2,
                                                                         alt1,
                                                                         exon_offset,
                                                                         0)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) #
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_alt1[event] = wiggle
        
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        alt1,
                                                                        alt2,
                                                                        exon_offset,
                                                                        0)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_alt2[event] = wiggle
                  
        three_upstream = pd.DataFrame(three_upstream).T
        five_alt1 = pd.DataFrame(five_alt1).T
        three_alt1 = pd.DataFrame(three_alt1).T
        five_alt2 = pd.DataFrame(five_alt2).T
        if combine_regions == False:
            return three_upstream, five_alt1, three_alt1, five_alt2
        else:
            ra = pd.concat([three_upstream,five_alt1,three_alt1,five_alt2],axis=1)
            ra.columns = range(0,ra.shape[1])
            return ra                           
                        
def create_se_matrix(annotation, density, exon_offset, intron_offset, is_scaled, combine_regions=True):
    print("creating se matrix for {} on {}".format(density.name,annotation))
    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}
    
            
    with open(annotation) as f:
        # f.next() # for title
        for line in f:
            if not line.startswith('#'):
                event = line.split('\t')[0]
                upstream_interval, interval, downstream_interval = Feature.SkippedExonFeature(event,'miso').get_bedtools()
                
                """ Deprecated
                upstream, se, downstream = event.split('@')
                
                upstream_interval = misc.create_bed_tool_from_miso_se(upstream)
                interval = misc.create_bed_tool_from_miso_se(se)
                downstream_interval = misc.create_bed_tool_from_miso_se(downstream)
                """
                """three prime upstream region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        interval,
                                                                        upstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                    
                wiggle = pd.Series(wiggle)
                # if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
        
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) 
                
                    # print("length of 3p upstream: {}".format(len(wiggle)))
                # else:
                #     wiggle = pd.Series([-1]*(intron_offset + exon_offset))
                three_upstream[event] = wiggle
                """five prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream_interval,
                                                                        interval,
                                                                        exon_offset,
                                                                        intron_offset)
                
                wiggle = pd.Series(wiggle)
                # if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle)
                
                    # print("length of 5p skipped: {}".format(len(wiggle)))
                # else:
                #     wiggle = pd.Series([-1]*(intron_offset + exon_offset))
                five_skipped[event] = wiggle
                """three prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                         downstream_interval,
                                                                         interval,
                                                                         exon_offset,
                                                                         intron_offset)
                wiggle = pd.Series(wiggle)
                # if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) #
                
                    # print("length of 3p skipped: {}".format(len(wiggle)))
                # else:
                #     wiggle = pd.Series([-1]*(intron_offset + exon_offset))
                three_skipped[event] = wiggle
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        interval,
                                                                        downstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                # if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                
                    # print("length of 5p downstream: {}".format(len(wiggle)))
                # else:
                #     wiggle = pd.Series([-1]*(intron_offset + exon_offset))
                five_downstream[event] = wiggle

        three_upstream = pd.DataFrame(three_upstream).T
        five_skipped = pd.DataFrame(five_skipped).T
        three_skipped = pd.DataFrame(three_skipped).T
        five_downstream = pd.DataFrame(five_downstream).T
        if combine_regions == False:
            return three_upstream, five_skipped, three_skipped, five_downstream
        else:
            ra = pd.concat([three_upstream,five_skipped,three_skipped,five_downstream],axis=1)
            ra.columns = range(0,ra.shape[1])
            # print("TYPE OF MATRIX: {}".format(type(ra)))
            return ra

# makes a raw count histogram of rbp for each position in each region
def make_hist_se(infile, outfile, hashing_val, l10p_cutoff, l2fc_cutoff, all_exons, exon_overhang, intron_overhang):
    try:
        region_types = ["upstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_skipped_exon",
                        "downstream_region_upstream_exon"]
        position_sum = {}
        count = 0
        with open(infile,'r') as f:
            for line in f:
                line = line.split('\t')
                chrom = line[0]
                pstart = int(line[1])
                pstop = int(line[2])
                l10p = float(line[3])
                l2fc = float(line[4])
                stra = line[5].strip()
                
                # correct bed files being 0-based, open ended
                pstart = pstart + 1
                
                if l10p < l10p_cutoff:
                    continue
                if l2fc < l2fc_cutoff:
                    continue
                
                x = int(pstart / hashing_val)
                y = int(pstop / hashing_val)

                # for each peak, find ALL regions that intersect it
                for region_type in region_types: # within a region
                    tmphash = {}
                    for i in range(x,y+1): # within a bin
                        for event in all_exons[chrom,stra,i,region_type]:
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            if pstop < int(exregstart): # pass if peak stop occurs before exon region start
                                continue
                            if pstart > int(exregstop): # pass if peak start occurs after exon region end
                                continue
                            tmphash[event] = 1 # otherwise peak falls within event region
                    for event in tmphash:
                        if stra == "+":
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart)) # peak start OR region start
                            end_val = min(int(pstop), int(exregstop)) # peak stop OR region stop
                            for j in range(start_val, end_val+1): # count intersecting positions between peak and region
                                relative_pos = j - int(exstart) # calculate relative position
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos) # count + 1 for the region
                        elif stra == '-':
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))
                            end_val = min(int(pstop), int(exregstop))
                            for j in range(start_val, end_val+1):
                                relative_pos = -1 * (j - int(exstart))
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos)
                        else:
                            print("strand error\n")
                    # we have a peak that maps to every region
                    
        # count from 0 to max
        current_pos = 0
        o = open(outfile,'w')
        for j in range(-exon_overhang, intron_overhang+1):
            if misc.exists(position_sum,"downstream_region_upstream_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_upstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang+1):
            if misc.exists(position_sum,"upstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-exon_overhang, intron_overhang+1):
            if misc.exists(position_sum,"downstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang+1):
            if misc.exists(position_sum,"upstream_region_downstream_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_downstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        o.close()
    except Exception as e:
        print(e)