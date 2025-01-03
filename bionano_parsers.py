from collections import defaultdict
from os.path import exists
import pandas as pd

#################################### parsing RefAligner Copy Number file ##########################################
def parse_rcmap(cmap_dir):
    cov = {}
    cov = defaultdict(lambda: {}, cov)
    cop = {}
    cop = defaultdict(lambda: {}, cop)
    with open(cmap_dir, 'r') as f:
        for line in f:
            if line.startswith("#"):
                head = line[1:].rstrip().rsplit()
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                chrom = fD['CMapId']
                pos = int(float(fD['Position']))
                cover = float(fD['Coverage'])
                copynumber = int(fD['CopyNumber'])
                cov[chrom][pos] = cover
                cop[chrom][pos] = copynumber
    return cov, cop # return two dictionary of dictionary which keys are chromosome number and keys are label position and return coverage/ copynumber
    #cov return coverage / #cop return copynumber

class BP: #Class of SV breakpoints 
    contig_id = '' # contig id the breakpoint is called
    direction1 = '' # Forward direction is + / Reverse direction is -
    direction2 = '' # Forward direction is + / Reverse direction is -
    pos1 = '' # Always lower genomic coordinate and lower chromosomes are sorted
    pos2 = ''
    chrom1 = ''
    chrom2 = ''
    line = ''
    type = ''

class SmapEntry: #Class of each entry in Smap file
    smap_id = ''
    q_id = ''
    ref_c_id1 = ''
    ref_c_id2 = ''
    ref_start = 0
    ref_end = 0
    query_start = 0
    query_end = 0
    confidence = 0
    xmap_id1 = ''
    xmap_id2 = ''
    sv_type = ''
    line = ''
    size = 0
    linkID = ''
    VAF = 0

#################################### parsing RefAligner smap file ##########################################
def parse_smap(smap_dir):
    with open(smap_dir, 'r') as f:
        bfb_count = {}
        bfb_count = defaultdict(lambda: [], bfb_count)
        breakpoints = []
        translocations = []
        for line in f:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                smap_entry = SmapEntry()
                smap_entry.ref_start = float(fD['RefStartPos'])
                smap_entry.ref_end = float(fD['RefEndPos'])
                smap_entry.xmap_id1 = int(fD['XmapID1'])
                smap_entry.xmap_id2 = int(fD['XmapID2'])
                smap_entry.q_id = int(fD['QryContigID'])
                smap_entry.ref_c_id1 = fD['RefcontigID1']
                smap_entry.ref_c_id2 = fD['RefcontigID2']
                smap_entry.smap_id = int(fD['SmapEntryID'])
                smap_entry.confidence = float(fD['Confidence'])
                smap_entry.query_start = float(fD['QryStartPos'])
                smap_entry.query_end = float(fD['QryEndPos'])
                smap_entry.sv_type = fD['Type']
                smap_entry.size =float(fD['SVsize'])
                smap_entry.line = line
                smap_entry.linkID = int(fD['LinkID'])
                # smap_entry.VAF = float(fD['VAF'])
                breakpoints.append(smap_entry)
    return breakpoints # return list of smap entries

class Segments: #Class represent each genomic segment
    id = ''
    chromosome = 0
    start = 0
    end = 0 
    width = 0
    type = ''
    fractional_cn = 0
    int_cn = 0
    conf = 0
    line = ''
    bp = [] # this list contains all SVs intersect middle of this genomic segment including start and end position of segment, we will use this then for spiliting this to new genomic segments

######################## parse CNV call ########################
def parse_cnvcall(cnvcall):
    segment_list = []
    all_seg = []
    with open(cnvcall, 'r') as f:
        for line in f:
            if line.startswith("#Id"):
                head = line.rstrip().rsplit()
            if not line.startswith('#'):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                segment = Segments()
                segment.id = int(fD['#Id'])
                segment.chromosome = fD['Chromosome']
                segment.start = float(fD['Start'])
                segment.end = float(fD['End'])
                if 'Width' in fD.keys():
                    segment.width = float(fD['Width'])
                else:
                    segment.width = float(fD['Size'])
                segment.type = fD['Type']
                segment.fractional_cn = float(fD['fractionalCopyNumber'])
                segment.int_cn = int(fD['CopyNumber'])
                segment.conf = float(fD['Confidence'])
                segment.line = line
                segment.bp = [segment.start, segment.end]
                if segment.width > 200000  and not segment.type.endswith('masked') and segment.conf>= 0.98: #Apply filters on CNV call masked region and segments legnth 200000bp is a limit of filtering segments
                    # if segment.width < 500000:
                    #     if len(segment_list) == 0:
                    #         segment_list.append(segment)
                    #     elif (segment.int_cn >= segment_list[-1].int_cn and segment.chromosome == segment_list[-1].chromosome) or segment.chromosome != segment_list[-1].chromosome: #If segment length is between 200 and 500 Kbp. if the CN is great than previouse segment length will add it( this prevent having small deletions between 200 and 500 Kbp)
                    #         segment_list.append(segment)
                    # else:# if segment length is greater than 500Kbp add it.
                        segment_list.append(segment)
                    # segment_list.append(segment)
                # elif segment.width > 200000  and segment.conf>= 0.99: #If segment is in the telomeric reason and we are high conf add it
                #     segment_list.append(segment)
                all_seg.append(segment)
    return segment_list, all_seg #return two lists of segment class. one the filtered one and one all of them. 

#################################### parsing Alignment xmap file ##########################################
def parse_xmap(xmapf):
    detailFields = ["XmapEntryID", "QryContigID", "RefContigID", "Orientation", "QryLen", "RefLen",
                    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "Alignment","Confidence"]
    numeric = ["QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos","Confidence"]
    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head, fields))
                alnstring = ")" + fD["Alignment"] + "("
                # xmapAln[fD["XmapEntryID"]] = alnstring.rsplit(")(")[1:-1]

                xmapPair[fD["XmapEntryID"]] = {x: fD[x] for x in detailFields}
                for keyword in numeric:
                    xmapPair[fD["XmapEntryID"]][keyword] = float(xmapPair[fD["XmapEntryID"]][keyword])

    return xmapPair #return  dict of dict which key is Xmapentry 

################################ Parse bed file ###########################
def parse_bed(bed_dir):
    l = []
    with open(bed_dir, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            chrom = int(line[0][3:])
            start = int(float(line[1]))
            end = int(float(line[2]))
            l.append([chrom, start, end])
    return l 
############################################# 
#This function parse the centromere region
def parse_centro(centro):
    r = {}
    r = defaultdict(lambda: [], r)
    with open(centro, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            key = line[0]
            pos1 = int(line[1])
            pos2 = int(line[2])
            r[key].append(pos1)
            r[key].append(pos2)
    return r


def parse_forbiden_region(forbiden):
    masked_regions = {}
    with open(forbiden, 'r') as file:
        lines = file.readlines()
    for line in lines[1:]:
        columns = line.strip().split('\t')
        if columns[0] == 'ChrX':
            columns[0] = 23
        elif columns[0] == 'ChrY':
            columns[0] = 24
        else:
            columns[0] = int(columns[0][3:])
        region = {
            'Chr': columns[0],
            'StartPos': int(columns[1]),
            'EndPos': int(columns[2]),
            'Type': columns[3]
        }
        if columns[0] not in masked_regions:
            masked_regions[columns[0]] = []
        masked_regions[columns[0]].append([int(columns[1]), int(columns[2]), columns[3]])
    return masked_regions


## ZJ
def smap_to_df(smap_filepath):
    with open(smap_filepath) as fp_read:
        line_count = 0
        for line in fp_read:
            line_count += 1
            if line.startswith('#h Smap'):
                break
        headers = line.replace('\n', '').split('\t')[1:]
    df = pd.read_csv(smap_filepath, sep='\t', skiprows=line_count+1, header=None, names=headers)
    return df

def cnv_to_df(cnv_filepath):
    skip = 2
    with open(cnv_filepath) as fp:
        line = fp.readline()
        if line.startswith('#Id'):
            skip = 0
    with open(cnv_filepath) as fp_read:
        for i in range(skip):
            fp_read.readline()
        headers = fp_read.readline().replace('\n', '').replace('#', '').split('\t')[1:]
    df = pd.read_csv(cnv_filepath, sep='\t', skiprows=skip+1, header=None, names=headers)
    return df

def edge_in_smap(smap_df, chrom1, chrom2, start, end, allowed_event_types, confidence_threshold, max_coord_distance):
    filtered_df = smap_df[(smap_df['RefcontigID1'] == chrom1) & (smap_df['RefcontigID2'] == chrom2)]
    if allowed_event_types:
        filtered_df = filtered_df[filtered_df['Type'].isin(allowed_event_types)]
    if confidence_threshold != -1:
        filtered_df = filtered_df[filtered_df['Confidence'] >= confidence_threshold]

    # ignore start/end orientation
    mask1 = abs(filtered_df['RefStartPos'] - start) + abs(filtered_df['RefEndPos'] - end) < max_coord_distance
    mask2 = abs(filtered_df['RefStartPos'] - end) + abs(filtered_df['RefEndPos'] - start) < max_coord_distance
    filtered_df = filtered_df[mask1 | mask2]
    return not filtered_df.empty

def region_reported_cn(cnv_df, chrom, start, end):
    """
    :param cnv_df:
    :param chrom: int (1-24)
    :param start:
    :param end:
    :return: average CN of the region, expected CN
    """
    ## assumes no overlap in row entry in the cnv_df
    filtered_df = cnv_df[cnv_df['Chromosome'] == chrom]
    filtered_df.sort_values(by='Start', ascending=True)

    ## figure out the expected count (mostly used for X)
    expected_cn = -1
    for index, row in filtered_df.iterrows():
        if 1 < row['fractionalCopyNumber'] < 2 and 'loss' in row['Type']:
            expected_cn = 2
            break
        elif 1 < row['fractionalCopyNumber'] < 2 and 'gain' in row['Type']:
            expected_cn = 1
            break
    if expected_cn == -1:
        print('assumed expected cn to be 2')
        expected_cn = 2

    ## tally average cn of the region from CNV file
    total_cnv_entries = []
    for index, row in filtered_df.iterrows():
        if start <= row['Start'] <= end <= row['End']:
            total_cnv_entries.append({'start': row['Start'], 'end': end, 'CN': row['fractionalCopyNumber']})
        elif row['Start'] <= start <= end <= row['End']:
            total_cnv_entries.append({'start': start, 'end': end, 'CN': row['fractionalCopyNumber']})
        elif row['Start'] <= start <= row['End'] <= end:
            total_cnv_entries.append({'start': start, 'end': row['End'], 'CN': row['fractionalCopyNumber']})
        elif start <= row['Start']<= row['End'] <= end:
            total_cnv_entries.append({'start': row['Start'], 'end': row['End'], 'CN': row['fractionalCopyNumber']})

    ## fill missing gap with no CNV
    total_cn = []
    c_start = start
    for cnv in total_cnv_entries:
        if cnv['start'] > c_start:
            total_cn.append({'start': c_start, 'end': cnv['start'] - 1, 'CN': expected_cn})
        total_cn.append(cnv)
        c_start = cnv['end'] + 1
    if not total_cnv_entries:
        total_cn.append({'start': start, 'end': end, 'CN': expected_cn})
    elif total_cnv_entries[-1]['end'] < end:
        total_cn.append({'start': total_cnv_entries[-1]['end'] + 1, 'end': end, 'CN': expected_cn})

    ## averaging CN
    cn_sum = 0
    for e in total_cn:
        cn_sum += e['CN'] * (e['end'] - e['start'] + 1)
    cn_sum /= (end - start + 1)

    return cn_sum, expected_cn
