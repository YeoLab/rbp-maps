'''
Created on May 19, 2016

@author: brianyee
'''

import sys
import os
import pandas as pd
import urllib, json, requests


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

def star_or_tophat(link):
    if os.path.basename(link).find('star') > -1:
        return 'STAR'
    if os.path.basename(link).find('tophat') > -1:
        return 'TOPHAT'
    return 'ambiguous'

def main(argv=None): # IGNORE:C0111
    
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    """
    infile: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    """
    parser.add_argument("-i", "--infile", dest="infile", help="input manifest file", required = True )
    parser.add_argument("-o", "--outfile", dest="outfile", help="outfile manifest file", required = True )
    parser.add_argument("-l", "--label", dest="label", help="lab name", required = False, default = 'encode-processing-pipeline')
    
    host = 'https://www.encodeproject.org'
    experiments = "https://www.encodeproject.org/experiments/"
    
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    labname = args.label
    
    samples = list()
    links = list()
    controls = list()
    error = list()
    
    X = pd.read_table(infile, index_col=0)
    
    for h in X['RNASEQ_ENCODEAccID'].dropna():
        url = experiments+h+"/?format=json"
        response = urllib.urlopen(url)
        data = json.loads(response.read())
        if 'code' in data.keys():
            next
        else:
            for i in range(0,len(data['files'])):
                # for key, value in data['files'][i].iteritems():
                #     print(key)
                # print(data['files'][i].iteritems())
                # print( data['files'][i]['replicate']['biological_replicate_number'])
                if ((host+data['files'][i]['href']).endswith('bam') & 
                (data['files'][i]['output_type'] == u'alignments') &
                (data['files'][i]['lab'][u'name'] == labname)):
                    if(list(X[X['RNASEQ_ENCODEAccID']==h]['RNASEQ_ControlENC'])[0] != data['possible_controls'][0]['accession']):
                        error.append("WARNING: CTRL VS AFF CHANGED: {}:{}:{}".format(
                            h,
                            list(X[X['RNASEQ_ENCODEAccID']==h]['RNASEQ_ControlENC'])[0],
                            data['possible_controls'][0]['accession']))
                    samples.append({'biological_replicate_number':data['files'][i]['replicate']['biological_replicate_number'],
                                    'submitted_expt_file_name':data['files'][i]['submitted_file_name'],
                                    'expt_file_name':os.path.basename(data['files'][i]['href']),
                                    'assembly':data['files'][i]['assembly'],
                                    'expt_url':host+data['files'][i]['href'],
                                    'RNASEQ_ControlENC':data['possible_controls'][0]['accession'],
                                    'RNASEQ_ENCODEAccID':h,
                                    'aligner': star_or_tophat(data['files'][i]['submitted_file_name'])
                                   })
                    links.append(host+data['files'][i]['href'])
            controls.append(data['possible_controls'][0]['accession'])
    samples_df = pd.DataFrame(samples)
    
    control_samples = list()
    error = list()
    for h in controls:
        url = experiments+h+"/?format=json"
        response = urllib.urlopen(url)
        data = json.loads(response.read())
        if 'code' in data.keys():
            next
        else:
            for i in range(0,len(data['files'])):
                # for key, value in data['files'][i].iteritems():
                #     print(key)
                # print(data['files'][i].iteritems())
                # print( data['files'][i]['replicate']['biological_replicate_number'])
                if ((host+data['files'][i]['href']).endswith('bam') & 
                (data['files'][i]['output_type'] == u'alignments') &
                (data['files'][i]['lab'][u'name'] == labname)):
                    control_samples.append({'biological_replicate_number':data['files'][i]['replicate']['biological_replicate_number'],
                              'submitted_ctrl_file_name':data['files'][i]['submitted_file_name'],
                              'ctrl_file_name':os.path.basename(data['files'][i]['href']),
                              'assembly':data['files'][i]['assembly'],
                              'ctrl_url':host+data['files'][i]['href'],
                              'aligner': star_or_tophat(data['files'][i]['submitted_file_name']),
                              'RNASEQ_ControlENC':h})
                    links.append(host+data['files'][i]['href'])
                    # controls[data['accession']] = data['possible_controls'][0]['accession']
    controls_df = pd.DataFrame(control_samples)
    
    merged_df = pd.merge(samples_df,controls_df,how='left',on=['RNASEQ_ControlENC','aligner','assembly','biological_replicate_number'])
    
    merged_df.to_csv(outfile,index=None,sep="\t")

if __name__ == '__main__':
    main()