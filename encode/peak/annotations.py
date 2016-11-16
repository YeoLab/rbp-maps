import collections
import pandas as pd
import sys

def read_alt_exons(alt_file):
    alt = open(alt_file)
    
def filter_summary(infile,included_outfile,excluded_outfile):
    oi = open(included_outfile,'w')
    oe = open(excluded_outfile,'w')
    with open(infile,'r') as f:
        for line in f:
            if line.startswith('event_name') == False:
                line = line.split('\t')
                event = line[0]
                inclusion = float(line[1])
                ci_low = float(line[2])
                ci_high = float(line[3])
                counts = line[6]
                sum_r1 = int(counts.split(',')[0].split(':')[1])
                try:
                    sum_r2 = int(counts.split(',')[1].split(':')[1])
                except IndexError:
                    sum_r2 = 0
                sums = sum_r1 + sum_r2
                if ((inclusion > 0.7) and (ci_high-ci_low < 0.3) and (sums > 10)):
                    oi.write("{0}\t{1}\n".format(event,"some_included_exon"))
                elif ((inclusion < 0.3) and (ci_high-ci_low < 0.3) and (sums > 10)):
                    oe.write("{0}\t{1}\n".format(event,"some_excluded_exon"))

def read_region_from_miso(event_line, event_type):
    upchr = ''
    upstart = ''
    upstop = ''
    upstr = ''
    if event_type == 'a3ss':
        event, ensg = event_line.split('\t')
        up, alt = event.split('@')
            
        upchr,upstart,upstop,upstr = up.split(':')
        alt = alt.split(':')
        exchr = alt[0]
        exstart, exstop = alt[1].split('|')
        
        exstr = alt[3]
        dnchr = exchr
        if(exstr == '+'):
            
            dnstart = int(exstop)
            exstop = int(exstop) - 1
            dnstop = alt[2]
        elif(exstr == '-'):
            
            dnstart = alt[2]
            dnstop = exstart
            exstart = int(exstart) + 1
        else:
            print("strand error")
        dnstr = exstr
    elif event_type == 'a5ss':
        # chr17:80009218:80008888|80009170:-@chr17:80008538:80008640:-
        # chr2:231331835:231331921|231331930:+@chr2:231333758:231333811:+
        event, ensg = event_line.split('\t')
        alt, dn = event.split('@')
        dnchr, dnstart, dnstop, dnstr = dn.split(':')
        alt = alt.split(':')
        exchr = alt[0]
        exstart, exstop = alt[2].split('|')
        exstr = alt[3]
        upchr = exchr
        if(exstr == '+'):
            upstart = alt[1]
            upstop = int(exstart) - 1
        elif(exstr == '-'):
            upstart = int(exstop) + 1
            upstop = alt[1]
        else:
            print("Strand error")
        upstr = exstr
            
        dnchr = exchr
    elif event_type == 'SE':
        event, ensg = event_line.split('\t')
        up, ex, dn = event.split('@')
            
        upchr,upstart,upstop,upstr = up.split(':')
        exchr,exstart,exstop,exstr = ex.split(':')
        dnchr,dnstart,dnstop,dnstr = dn.split(':')
        
    upstream_exon = [upchr, upstart, upstop, upstr]
    excluded_exon = [exchr, exstart, exstop, exstr]
    downstream_exon = [dnchr, dnstart, dnstop, dnstr]
    
    return upstream_exon, excluded_exon, downstream_exon

def read_four_region_miso(infile, hashing_val, event_type, exon_overhang, intron_overhang):
    all_exons = collections.defaultdict(list)
    with open(infile,'r') as f:
        next(f) # header
        for line in f:
            upstream, excluded, downstream = read_region_from_miso(line, event_type)
            upchr, upstart, upstop, upstr = upstream
            exchr, exstart, exstop, exstr = excluded
            dnchr, dnstart, dnstop, dnstr = downstream
            """print("up: {}, {}, {}, {}".format(upchr, upstart, upstop, upstr))
            print("ex: {}, {}, {}, {}".format(exchr, exstart, exstop, exstr))
            print("dn: {}, {}, {}, {}".format(dnchr, dnstart, dnstop, dnstr))
            
            return 0"""
            if (exstr == '+'):
                # upstream region of skipped exons
                if event_type == 'a5ss':
                    upstream_region_start = int(upstop) 
                else:
                    upstream_region_start = max(int(upstop)+1,int(exstart)-intron_overhang) # (end of upstream exon) OR (start - 500)
                
                upstream_region_stop = min(int(exstop), int(exstart)+exon_overhang) # (end of skipped exon) OR (start of skipped exon + 50)
                
                ex_upstreamregion = "{0}:{1}:{2}:{3}:{4}".format(exchr,
                                                                 upstream_region_start,
                                                                 exstart,
                                                                 upstream_region_stop,
                                                                 exstr)
                
                x = int(upstream_region_start / hashing_val)
                y = int(upstream_region_stop / hashing_val)
                
                for i in range(x,y+1):
                    all_exons[exchr,exstr,i,'upstream_region_skipped_exon'].append(ex_upstreamregion)
                    #print(all_exons)
                
                # upstream region of downstream exon
                dnex_region_start = max(int(exstop)+1,int(dnstart)-intron_overhang) # end of skipped exon OR start of downstream exon - 500
                dnex_region_stop = min(int(dnstop),int(dnstart)+exon_overhang) # end of downstream exon OR start of downstream exon + 50
                dnex_upstreamregion = "{0}:{1}:{2}:{3}:{4}".format(dnchr,
                                                                   dnex_region_start,
                                                                   dnstart,
                                                                   dnex_region_stop,
                                                                   dnstr)
                x = int(dnex_region_start / hashing_val)
                y = int(dnex_region_stop / hashing_val)
                
                for i in range(x,y+1):
                    all_exons[exchr,exstr,i,'upstream_region_downstream_exon'].append(dnex_upstreamregion)
                    
                # downstream region of skipped exons
                downstream_region_start = max(int(exstart),int(exstop)-exon_overhang) # beginning of skipped exon OR end of skipped exon - 50
                downstream_region_stop = min(int(dnstart)-1,int(exstop)+intron_overhang) # start of downstream exon OR end of skipped exon + 500
                ex_downstreamregion = "{0}:{1}:{2}:{3}:{4}".format(exchr,
                                                                   downstream_region_start,
                                                                   exstop,
                                                                   downstream_region_stop,
                                                                   exstr)
                x = int(downstream_region_start / hashing_val)
                y = int(downstream_region_stop / hashing_val)
                
                for i in range(x,y+1):
                    all_exons[exchr,exstr,i,'downstream_region_skipped_exon'].append(ex_downstreamregion)
                    
                # downstream region of upstream exon
                upex_region_start = max(int(upstart),int(upstop)-exon_overhang) # start of upstream exon OR end of upstream exon - 50
                upex_region_stop = min(int(exstart)-1,int(upstop)+intron_overhang) # start of skipped exon OR end of upstream exon + 500
                upex_downstreamregion = "{0}:{1}:{2}:{3}:{4}".format(upchr,
                                                                     upex_region_start,
                                                                     upstop,
                                                                     upex_region_stop,
                                                                     upstr)
                
                x = int(upex_region_start / hashing_val)
                y = int(upex_region_stop / hashing_val)
                
                for i in range(x,y+1):
                    all_exons[exchr,exstr,i,'downstream_region_upstream_exon'].append(upex_downstreamregion)
            elif (exstr == '-'):
                # calculate upstream_region_skipped_exon, or UPSTREAM region of SKIPPED EXON
                if(event_type == 'a5ss'):
                    upstream_region_start = int(exstop)
                else:
                    upstream_region_start = min(int(upstart)-1,int(exstop)+intron_overhang)
                upstream_region_stop = max(int(exstart),int(exstop)-exon_overhang)
                
                ## upstream_region_stop and upstream_region_start are flipped due to strand
                ex_upstreamregion = "{0}:{1}:{2}:{3}:{4}".format(exchr,
                                                                 upstream_region_stop,
                                                                 exstop,
                                                                 upstream_region_start,
                                                                 exstr)
                x = int(upstream_region_stop / hashing_val)
                y = int(upstream_region_start / hashing_val)
                
                for i in range(x,y+1):
                    all_exons[exchr,exstr,i,'upstream_region_skipped_exon'].append(ex_upstreamregion)
                
                # upstream region of downstream exon
                dnex_region_start = min(int(exstart)-1,int(dnstop)+intron_overhang)
                dnex_region_stop = max(int(dnstart),int(dnstop)-exon_overhang)
                ## dnex_region_start and dnex_region_stop are flipped due to strand
                dnex_upstreamregion = "{0}:{1}:{2}:{3}:{4}".format(dnchr,
                                                                   dnex_region_stop,
                                                                   dnstop,
                                                                   dnex_region_start,
                                                                   dnstr)
                x = int(dnex_region_stop / hashing_val)
                y = int(dnex_region_start / hashing_val)

                for i in range(x,y+1):
                    all_exons[exchr,exstr,i,'upstream_region_downstream_exon'].append(dnex_upstreamregion)
                
                # downstream region of skipped exon 
                downstream_region_start = min(int(exstart)+exon_overhang,int(exstop))
                downstream_region_stop = max(int(exstart)-intron_overhang,int(dnstop)+1)
                ex_downstreamregion = "{0}:{1}:{2}:{3}:{4}".format(exchr,
                                                                   downstream_region_stop,
                                                                   exstart,
                                                                   downstream_region_start,
                                                                   exstr)
                x = int(downstream_region_stop / hashing_val)
                y = int(downstream_region_start / hashing_val)
                
                for i in range(x,y+1):
                    all_exons[exchr,exstr,i,'downstream_region_skipped_exon'].append(ex_downstreamregion)
                
                # downstream region of upstream exon    
                upex_region_start = min(int(upstart)+exon_overhang,int(upstop))
                upex_region_stop = max(int(upstart)-intron_overhang,int(exstop)+1)
                upex_downstreamregion = "{0}:{1}:{2}:{3}:{4}".format(upchr,
                                                                     upex_region_stop,
                                                                     upstart,
                                                                     upex_region_start,
                                                                     upstr)
                x = int(upex_region_stop / hashing_val)
                y = int(upex_region_start / hashing_val)
                
                for i in range(x,y+1):
                    all_exons[exchr,exstr,i,'downstream_region_upstream_exon'].append(upex_downstreamregion)
            else:
                print("strand error")
            # print all_exons
            # return 0
    return all_exons
def main():
    filter_summary("SE.miso_summary",
                   "SE_included.filtered.miso_summary",
                   "SE_excluded.filtered.miso_summary")
if __name__ == '__main__':
    sys.exit(main())