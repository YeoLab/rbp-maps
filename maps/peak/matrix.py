import misc


# makes a raw count histogram of rbp for each position in each region
def make_hist_se(infile, outfile, hashing_val, l10p_cutoff, l2fc_cutoff, all_exons, exon_overhang, intron_overhang):
    try:
        region_types = ["upstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_skipped_exon",
                        "downstream_region_upstream_exon"]
        position_sum = {}
        count = []
        with open(infile, 'r') as f:
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
                for region_type in region_types:  # within a region
                    tmphash = {}
                    for i in range(x, y + 1):  # within a bin
                        for event in all_exons[chrom, stra, i, region_type]:
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            if pstop < int(exregstart):  # pass if peak stop occurs before exon region start
                                continue
                            if pstart > int(exregstop):  # pass if peak start occurs after exon region end
                                continue
                            tmphash[event] = 1  # otherwise peak falls within event region
                    for event in tmphash:
                        if stra == "+":
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))  # peak start OR region start
                            end_val = min(int(pstop), int(exregstop))  # peak stop OR region stop
                            for j in range(start_val,
                                           end_val + 1):  # count intersecting positions between peak and region
                                relative_pos = j - int(exstart)  # calculate relative position
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                                   region_type,
                                                                                   relative_pos)  # count + 1 for the region
                        elif stra == '-':
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))
                            end_val = min(int(pstop), int(exregstop))
                            for j in range(start_val, end_val + 1):
                                relative_pos = -1 * (j - int(exstart))
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                                   region_type,
                                                                                   relative_pos)
                        else:
                            print("strand error\n")

        # count from 0 to max
        current_pos = 0
        o = open(outfile, 'w')
        for j in range(-exon_overhang, intron_overhang + 1):
            if misc.exists(position_sum, "downstream_region_upstream_exon", j):
                o.write("{}\n".format(position_sum["downstream_region_upstream_exon", j]))
                count.append(position_sum["downstream_region_upstream_exon", j])
            else:
                o.write("{}\n".format(0))
                count.append(0)
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang + 1):
            if misc.exists(position_sum, "upstream_region_skipped_exon", j):
                o.write("{}\n".format(position_sum["upstream_region_skipped_exon", j]))
                count.append(position_sum["upstream_region_skipped_exon", j])
            else:
                o.write("{}\n".format(0))
                count.append(0)
            current_pos = current_pos + 1
        for j in range(-exon_overhang, intron_overhang + 1):
            if misc.exists(position_sum, "downstream_region_skipped_exon", j):
                o.write("{}\n".format(position_sum["downstream_region_skipped_exon", j]))
                count.append(position_sum["downstream_region_skipped_exon", j])
            else:
                o.write("{}\n".format(0))
                count.append(0)
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang + 1):
            if misc.exists(position_sum, "upstream_region_downstream_exon", j):
                o.write("{}\n".format(position_sum["upstream_region_downstream_exon", j]))
                count.append(position_sum["upstream_region_downstream_exon", j])
            else:
                o.write("{}\n".format(0))
                count.append(0)
            current_pos = current_pos + 1
        o.close()
        return count
    except Exception as e:
        print(e)


# makes a raw count histogram of rbp for each position in each region
def make_hist_a3ss(infile, outfile, hashing_val, l10p_cutoff, l2fc_cutoff, all_exons, exon_overhang, intron_overhang):
    try:
        region_types = ["upstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_skipped_exon",
                        "downstream_region_upstream_exon"]
        position_sum = {}
        count = 0
        with open(infile, 'r') as f:
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
                for region_type in region_types:  # within a region
                    tmphash = {}
                    for i in range(x, y + 1):  # within a bin
                        for event in all_exons[chrom, stra, i, region_type]:
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            if pstop < int(exregstart):  # pass if peak stop occurs before exon region start
                                continue
                            if pstart > int(exregstop):  # pass if peak start occurs after exon region end
                                continue
                            tmphash[event] = 1  # otherwise peak falls within event region
                    for event in tmphash:
                        if stra == "+":
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))  # peak start OR region start
                            end_val = min(int(pstop), int(exregstop))  # peak stop OR region stop
                            for j in range(start_val,
                                           end_val + 1):  # count intersecting positions between peak and region
                                relative_pos = j - int(exstart)  # calculate relative position
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                                   region_type,
                                                                                   relative_pos)  # count + 1 for the region
                        elif stra == '-':
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))
                            end_val = min(int(pstop), int(exregstop))
                            for j in range(start_val, end_val + 1):
                                relative_pos = -1 * (j - int(exstart))
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                                   region_type,
                                                                                   relative_pos)
                        else:
                            print("strand error\n")

        # count from 0 to max
        current_pos = 0
        o = open(outfile, 'w')
        for j in range(-exon_overhang, intron_overhang + 1):
            if misc.exists(position_sum, "downstream_region_upstream_exon", j):
                o.write("{}\n".format(position_sum["downstream_region_upstream_exon", j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang + 1):
            if misc.exists(position_sum, "upstream_region_skipped_exon", j):
                o.write("{}\n".format(position_sum["upstream_region_skipped_exon", j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-exon_overhang, 1):
            if misc.exists(position_sum, "downstream_region_skipped_exon", j):
                o.write("{}\n".format(position_sum["downstream_region_skipped_exon", j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(0, exon_overhang + 1):
            if misc.exists(position_sum, "upstream_region_downstream_exon", j):
                o.write("{}\n".format(position_sum["upstream_region_downstream_exon", j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        o.close()
    except Exception as e:
        print(e)


def make_hist_a5ss(infile, outfile, hashing_val, l10p_cutoff, l2fc_cutoff, all_exons, exon_overhang, intron_overhang):
    try:
        # region_types = ["ex_up","dnex","ex_dn","upex"] # former annotation
        region_types = ["upstream_region_skipped_exon",
                        "downstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_upstream_exon"]
        position_sum = {}
        count = 0
        with open(infile, 'r') as f:
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
                for region_type in region_types:  # within a region
                    tmphash = {}
                    for i in range(x, y + 1):  # within a bin
                        for event in all_exons[chrom, stra, i, region_type]:
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            if pstop < int(exregstart):  # pass if peak stop occurs before exon region start
                                continue
                            if pstart > int(exregstop):  # pass if peak start occurs after exon region end
                                continue
                            tmphash[event] = 1  # otherwise peak falls within event region
                    for event in tmphash:
                        if stra == "+":
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))  # peak start OR region start
                            end_val = min(int(pstop), int(exregstop))  # peak stop OR region stop
                            for j in range(start_val,
                                           end_val + 1):  # count intersecting positions between peak and region
                                relative_pos = j - int(exstart)  # calculate relative position
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                                   region_type,
                                                                                   relative_pos)  # count + 1 for the region
                        elif stra == '-':
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))
                            end_val = min(int(pstop), int(exregstop))
                            for j in range(start_val, end_val + 1):
                                relative_pos = -1 * (j - int(exstart))
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                                   region_type,
                                                                                   relative_pos)
                        else:
                            print("strand error\n")

        # count from 0 to max
        current_pos = 0
        o = open(outfile, 'w')
        for j in range(-exon_overhang, 1):
            if misc.exists(position_sum, "downstream_region_upstream_exon", j):
                o.write("{}\n".format(position_sum["downstream_region_upstream_exon", j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        # o.write("STOP")
        for j in range(0, exon_overhang + 1):
            if misc.exists(position_sum, "upstream_region_skipped_exon", j):
                o.write("{}\n".format(position_sum["upstream_region_skipped_exon", j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        # o.write("STOP")
        for j in range(-exon_overhang, intron_overhang + 1):
            if misc.exists(position_sum, "downstream_region_skipped_exon", j):
                o.write("{}\n".format(position_sum["downstream_region_skipped_exon", j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang + 1):
            if misc.exists(position_sum, "upstream_region_downstream_exon", j):
                o.write("{}\n".format(position_sum["upstream_region_downstream_exon", j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        o.close()
    except Exception as e:
        print(e)
