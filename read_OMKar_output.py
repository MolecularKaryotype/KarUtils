from .forbidden_region_processing import *
from .utils import *
from collections import defaultdict
import pandas as pd

def read_OMKar_output_to_path(OMKar_output_file, forbidden_region_file):
    path_list, index_dict = read_OMKar_output(OMKar_output_file, return_segment_dict=True)
    label_path_with_forbidden_regions(path_list, forbidden_region_file)
    rotate_and_bin_path(path_list, forbidden_region_file)
    report_centromere_anomaly(path_list)
    return index_dict, path_list

def group_segments_by_chr(segment_dict):

    # Initialize defaultdict to hold lists
    grouped_by_chr = defaultdict(list)

    # Iterate over the people and group by age
    for index,segment in segment_dict.items():
        chr = segment.chr_name
        grouped_by_chr[chr].append(segment)

    # Convert defaultdict back to a regular dictionary if needed
    grouped_by_chr = dict(grouped_by_chr)
    return grouped_by_chr

def all_segments_continuous(segments):
    # Iterate through the list of segments
    for i in range(len(segments) - 1):
        if not segments[i].continuous(segments[i + 1]):
            return False  
    return True 

def check_continous(segment_dict):
    groups = group_segments_by_chr(segment_dict)
    for chr, group in groups.items():
        if all_segments_continuous(group) == False:
            return False
    return True

# Validate the omkar output is correct
def validate_OMKar_output(path_list, segment_dict):
    checker = True
    
    for index,segment in segment_dict.items():
        if not segment.direction():
            checker = False
    # Implement the check_continous function to group by chr and check each segment.
    checker = check_continous(segment_dict=segment_dict)
                           
    return checker

def concatenate_segments(segments):
    return ' '.join(segment.kt_index for segment in segments)


def post_process_function(path_list, segment_list):
    tmp_path_list = []
    for path in path_list:
        tmp_path = path.duplicate()
        tmp_path_list.append(tmp_path)
    # label_path_with_forbidden_regions(tmp_path_list, forbidden_region_file)
    rotated_path_idx = rotate_and_bin_path(tmp_path_list, return_rotated_idx=True)
    for path_idx in rotated_path_idx:
        path = path_list[path_idx]
        path.reverse()

    merged_segments = []
    path_segment_dict = {}
    
    for path in path_list:
        new_concate_segments = concatenate_segments(path.linear_path.segments)
        path_segment_dict[path.path_name] = new_concate_segments
    duplicated_v = repeated_segments(path_segment_dict)
    
    #For each path that has segments to be merged, update segment list
    for segments, paths in duplicated_v.items():
        for path_name in paths:
            path = find_path(path_list, path_name)
            path.linear_path.segments = merge_segments(path.linear_path.segments)
    segments = []
    for path in path_list:
        segments.extend(path.linear_path.segments)
    sorted_segments = list(set(segments))
    sorted_segments = sorted(sorted_segments, key=sort_segment)
    for index, seg in enumerate(sorted_segments):
        seg_dir = seg.kt_index[-1]
        seg.kt_index = str(index)+seg_dir

def sort_segment(segment):
    return int(segment.kt_index[:-1])

def find_path(path_list, path_name):
    for path in path_list:
        if path.path_name == path_name:
            return path
    return None

def merge_segments(segment_list):
    #TODO Make sure that the list of segments are valid
    segment_direction = segment_list[0].kt_index[-1]
    segment = Segment(segment_list[0].chr_name, segment_list[0].start, segment_list[-1].end, "OMKar_unlabeled")
    segment.kt_index = segment_list[0].kt_index
    return [segment]

def repeated_segments(path_segment_dict):
    # Reverse mapping dictionary
    reverse_dict = {}

    # Populate reverse_dict with segments as keys and paths as values
    for path, segments in path_segment_dict.items():
        if segments not in reverse_dict:
            reverse_dict[segments] = []
        reverse_dict[segments].append(path)

    # Find duplicates
    duplicates = {seg: paths for seg, paths in reverse_dict.items() if len(paths) >1 and len(seg.split(" ")) > 1}
    return duplicates


    
    # merge_identical_paths(path_list)


def write_OMKar_to_file(path_list, segment_list, filename):
    return None

def read_OMKar_output(file, return_segment_dict=False):
    segment_dict = {}
    path_list = []
    with open(file) as fp_read:
        fp_read.readline()

        for line in fp_read:
            line = line.replace("\n", "").split('\t')
            # documenting segments
            if line[0] == "Segment":
                chr_name = str(line[2])
                if chr_name == "23":
                    chr_name = "X"
                elif chr_name == "24":
                    chr_name = "Y"
                chr_name = "Chr" + chr_name
                # start = int(line[3].split(".")[0])
                # end = int(line[4].split(".")[0])
                start = round_half_up(float(line[3]))
                end = round_half_up(float(line[4]))
                segment_dict[int(line[1])] = Segment(chr_name, start, end, "OMKar_unlabeled")
            elif line[0].startswith("Path"):
                # print(line)
                line = line[0].split(" = ")
                path_name = line[0]
                line = line[1]
                line = line.split(" ")
                path_segments = []
                for segment_index_itr in line:
                    if len(segment_index_itr) == 0:
                        break
                    direction = segment_index_itr[-1]
                    segment_index_itr = int(segment_index_itr[:-1])
                    new_segment = segment_dict[segment_index_itr].duplicate()
                    new_segment.kt_index = str(segment_index_itr)
                    new_segment.kt_index += '+'
                    if direction == "+":
                        path_segments.append(new_segment)
                    elif direction == "-":
                        new_segment.invert()
                        path_segments.append(new_segment)
                    else:
                        # print(direction)
                        raise ValueError("direction must be + or -")
                path_list.append(Path(Arm(path_segments, "solved_path"), path_name))

    if return_segment_dict:
        return path_list, segment_dict
    else:
        return path_list


def read_OMKar_to_indexed_list(OMKar_output_file, forbidden_region_file='Metadata/acrocentric_telo_cen.bed'):
    path_list, index_dict = read_OMKar_output(OMKar_output_file, return_segment_dict=True)
    ## extract which path to rotate and rotate, without splitting segments
    tmp_path_list = []
    for path in path_list:
        tmp_path = path.duplicate()
        tmp_path_list.append(tmp_path)
    label_path_with_forbidden_regions(tmp_path_list, forbidden_region_file)
    rotated_path_idx = rotate_and_bin_path(tmp_path_list, forbidden_region_file, return_rotated_idx=True)

    for path_idx, path in enumerate(path_list):
        if path_idx in rotated_path_idx:
            rotate_path(path)

    ## extract path's characterized chr
    path_chrs = []
    for path in tmp_path_list:
        path_chrs.append(path.path_chr)

    ## match and translate to indexing
    segment_dict = reverse_dict(index_dict)
    indexed_lists = []

    for path in path_list:
        indexed_list = []
        segments = path.linear_path.segments
        for segment in segments:
            if segment in segment_dict:
                indexed_list.append(str(segment_dict[segment]) + "+")
            else:
                segment_copy = segment.duplicate()
                segment_copy.invert()
                if segment_copy in segment_dict:
                    indexed_list.append(str(segment_dict[segment_copy]) + "-")
                else:
                    raise RuntimeError('segment_dict not complete')
        indexed_lists.append(indexed_list)

    segment_size_dict = {}
    for typed_seg, index_seg in segment_dict.items():
        segment_size_dict[str(index_seg)] = len(typed_seg)
    return indexed_lists, path_chrs, segment_dict, segment_size_dict


def generate_wt_from_OMKar_output(segment_to_index_dict):
    sorted_segments = sorted(list(segment_to_index_dict.keys()))
    wt_indexed_paths = {}
    c_chr = 'Chr1'
    c_path = []
    for segment in sorted_segments:
        seg_chr = segment.chr_name
        seg_index = segment_to_index_dict[segment]
        if seg_chr != c_chr:
            # new chr section entered
            wt_indexed_paths[c_chr] = c_path
            c_chr = seg_chr
            c_path = []
        c_path.append(str(seg_index) + '+')
    wt_indexed_paths[c_chr] = c_path

    return wt_indexed_paths


def rotate_and_bin_path(path_list, forbidden_region_file='Metadata/acrocentric_telo_cen.bed', return_rotated_idx=False):
    """
    only works if each path contains exactly one centromere, OW will bin according to t1+t2+centromere percentage,
    if still no, will bin according to overall chr-content percentage
    will mark path accordingly if centromere anomaly exists
    :param forbidden_region_file:
    :param path_list:
    :return: path_list
    """
    rotated_path_idx = []
    # isolate centromere
    # forbidden_region_path = Path(read_forbidden_regions(forbidden_region_file), 'forbidden_regions', 'forbidden_regions')
    label_path_with_forbidden_regions(path_list, forbidden_region_file)
    # for path in path_list:
    #     path.generate_mutual_breakpoints(other_path=forbidden_region_path, mutual=False)

    # get centromere, rotate if backward, and bin path
    for path_idx, path in enumerate(path_list):
        # print(path.path_name)
        path_centromere = []
        path_telomeres = []
        for segment_itr in path.linear_path.segments:
            if 'centromere' in segment_itr.segment_type:
                path_centromere.append(segment_itr.duplicate())
            elif 'telomere' in segment_itr.segment_type:
                path_telomeres.append(segment_itr.duplicate())

        path_centromere_arm = Arm(path_centromere, 'centromeres')
        path_centromere_arm.merge_breakpoints()

        if len(path_centromere_arm.segments) >= 1:
            centromere_set = set()
            for cen_seg in path_centromere_arm.segments:
                centromere_set.add(cen_seg.chr_name)
            if len(centromere_set) > 1:
                path.path_chr = '-multiple centromeres ({}), highest representation: {}'.format(centromere_set, get_highest_represented_chr(path_centromere))
            else:
                path.path_chr = path_centromere_arm.segments[0].chr_name
            if not highest_represented_direction(path_centromere):
                rotated_path_idx.append(path_idx)
                rotate_path(path)
        else:
            # no centromere segment detected, assume only q arm remains
            # take the first and the last segment, rotate chr if last segment index > first AND they are on the same chr
            path.path_chr = "-no centromere, highest representation: " + get_highest_represented_chr(path.linear_path.segments)
            first_segment = path.linear_path.segments[0]
            last_segment = path.linear_path.segments[-1]

            if first_segment.chr_name != last_segment.chr_name:
                # search for the last segment that is the same chr origin as the first segment (cont.)
                next_idx = 0
                while True:
                    if next_idx + 1 < len(path.linear_path.segments) and path.linear_path.segments[next_idx + 1].chr_name == first_segment.chr_name:
                        next_idx += 1
                    else:
                        break
                last_segment_same_chr = path.linear_path.segments[next_idx]
                forward_delta = last_segment_same_chr.end - first_segment.start

                # search for the first segment that is the same chr origin as the last segment (cont.)
                previous_idx = len(path.linear_path.segments) - 1
                while True:
                    if previous_idx - 1 >= 0 and path.linear_path.segments[previous_idx - 1].chr_name == last_segment.chr_name:
                        previous_idx -= 1
                    else:
                        break
                first_segment_same_chr = path.linear_path.segments[previous_idx]
                reverse_delta = last_segment.end - first_segment_same_chr.start

                # take major representation
                if forward_delta + reverse_delta < 0:
                    rotated_path_idx.append(path_idx)
                    rotate_path(path)
            else:
                if first_segment.start > last_segment.end:
                    rotated_path_idx.append(path_idx)
                    rotate_path(path)
    if return_rotated_idx:
        return rotated_path_idx


def count_chr_number(binned_path_list):
    chr_count = {}
    for path in binned_path_list:
        path_chr = path.path_chr
        if ": " in path_chr:
            path_chr = path_chr.split(': ')[-1]

        if path_chr not in chr_count:
            chr_count[path_chr] = 1
        else:
            chr_count[path_chr] += 1

    for i in range(1, 23):
        if f"Chr{i}" in chr_count and chr_count[f"Chr{i}"] == 2:
            chr_count.pop(f"Chr{i}")
        elif f"Chr{i}" not in chr_count:
            chr_count[f"Chr{i}"] = 0
    return chr_count


def report_centromere_anomaly(path_list):
    for path in path_list:
        if "no centromere" in path.path_chr or "multiple centromere" in path.path_chr:
            print(path.get_path_notes())


def get_highest_represented_chr(segment_list):
    tally = {}
    for segment in segment_list:
        if segment.chr_name in tally:
            tally[segment.chr_name] += len(segment)
        else:
            tally[segment.chr_name] = len(segment)

    max_count = -1
    max_count_chr = None
    for key in tally:
        if tally[key] > max_count:
            max_count = tally[key]
            max_count_chr = key
    return max_count_chr


def highest_represented_direction(segment_list):
    forward_len = 0
    backward_len = 0
    for segment in segment_list:
        if segment.direction():
            forward_len += len(segment)
        else:
            backward_len += len(segment)
    if forward_len >= backward_len:
        return True
    else:
        return False


def rotate_path(input_path):
    segment_list = input_path.linear_path.segments
    segment_list.reverse()
    for segment_itr in segment_list:
        segment_itr.invert()
    input_path.linear_path.segments = segment_list

# def cmd_centromere_anomaly():
#     # TODO: verify errors before restoring
#     import argparse
#     parser = argparse.ArgumentParser(description="dicentromeric and acentromeric checker")
#     parser.add_argument("--file", type=str, dest='omkar_file', help="file path to OMKar's solved path")
#     args = parser.parse_args()
#
#     path_list = read_OMKar_output(args.omkar_file)
#     for path in path_list:
#         path.linear_path.merge_breakpoints()
#     path_list = rotate_and_bin_path(path_list, "Metadata/merged_forbidden_regions_unique.bed")
#     report_centromere_anomaly(path_list)


### pending updates
def get_segments_by_type(forbidden_region_file, segment_type):
    masking_arm = read_forbidden_regions(forbidden_region_file)
    centromere_segments = []
    for segment in masking_arm.segments:
        if segment.segment_type == segment_type:
            centromere_segments.append(segment)
    return Arm(centromere_segments, segment_type)


def bin_path_by_chr_content(input_path):
    tally = {}
    for segment in input_path.linear_path.segments:
        if segment.chr_name in tally:
            tally[segment.chr_name] += len(segment)
        else:
            tally[segment.chr_name] = len(segment)

    max_count = -1
    max_count_chr = None
    for key in tally:
        if tally[key] > max_count:
            max_count = tally[key]
            max_count_chr = key
    return max_count_chr


def read_bed_file(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df


def bed_similar_sv_edge(bed_df, chromosome, pos1, pos2, approx_distance, reverse_search=True):
    chr_mask = bed_df['chromosome'] == chromosome
    if reverse_search:
        mask1 = (abs(bed_df['start'] - pos1) + abs(bed_df['end'] - pos2)) < approx_distance
        mask2 = (abs(bed_df['start'] - pos2) + abs(bed_df['end'] - pos1)) < approx_distance
        return bed_df[chr_mask & (mask1 | mask2)]
    else:
        mask1 = (abs(bed_df['start'] - pos1) + abs(bed_df['end'] - pos2)) < approx_distance
        return bed_df[chr_mask & mask1]


def test():
    path_list = read_OMKar_output("/media/zhaoyang-new/workspace/KarSim/KarComparator/new_data_files/OMKar_testbuild3/23X_15q26_overgrowth_r1.1.txt")
    rotate_and_bin_path(path_list, "Metadata/merged_forbidden_regions_unique.bed")
    # report_centromere_anomaly(path_list)
    for path in path_list:
        print(path)


def test_read_OMKar_output():
    path_list, segment_list = read_OMKar_output("sample_input/23Y_Cri_du_Chat_r1.1.txt", return_segment_dict=True)
    if validate_OMKar_output(path_list,segment_dict=segment_list) == True:
        post_process_function(path_list,segment_list)

def test_read_OMKar_to_path():
    idx_dict, path_list = read_OMKar_output_to_path("sample_input/23Y_Cri_du_Chat_r1.1.txt", "Metadata/acrocentric_telo_cen.bed")
    for path in path_list:
        print(path)


def test_output_index_list():
    # idx_dict, path_list = read_OMKar_output_to_path("sample_input/23Y_Cri_du_Chat_r1.1.txt", "Metadata/acrocentric_telo_cen.bed")
    indexed_lists, segment_dict, _ = read_OMKar_to_indexed_list("sample_input/23Y_Cri_du_Chat_r1.1.txt", "Metadata/acrocentric_telo_cen.bed")
    for lst in indexed_lists:
        print(lst)
    print('wt_list')
    wt_lists = generate_wt_from_OMKar_output(segment_dict)
    for lst in wt_lists:
        print(lst)


if __name__ == "__main__":
    # read_bed_file('/media/zhaoyang-new/workspace/keyhole/0717_output/510/510_SV.bed')
    test_read_OMKar_output()
    
