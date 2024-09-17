from cgitb import small

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

#################################MK-VALIDATION################################

def group_segments_by_chr(segment_dict):
    grouped_by_chr = defaultdict(list)

    # make sure the segment are enuemrated in order
    sorted_idx = sorted(segment_dict.keys(), key=int)
    for idx in sorted_idx:
        segment = segment_dict[idx]
        grouped_by_chr[segment.chr_name].append(segment)

    grouped_by_chr = dict(grouped_by_chr)
    return grouped_by_chr

def all_segments_continuous(segments, allowance):
    for i in range(len(segments) - 1):
        if not segments[i].validate_continuity(segments[i + 1], allowance):
            print(f"not continuous: {segments[i]}, {segments[i + 1]}")
            return False
    return True

def check_continous(segment_dict, allowance=5):
    groups = group_segments_by_chr(segment_dict)
    for chrom, segments in groups.items():
        if not all_segments_continuous(segments, allowance):
            return False
    return True

def check_spanning(segment_dict, forbidden_region_file, allowance=5):
    groups = group_segments_by_chr(segment_dict)
    nonforbidden_boundaries = get_prefix_suffix_forbidden_boundaries(forbidden_region_file=forbidden_region_file)
    for chrom, segments in groups.items():
        chrom_nonforbidden_boundaries = nonforbidden_boundaries[chrom]
        ordered_segments = [seg.duplicate() for seg in segments]
        for seg in ordered_segments:
            if not seg.direction():
                seg.invert()
        current_pos = chrom_nonforbidden_boundaries['start']
        for seg in ordered_segments:
            if seg.start - current_pos >= allowance:
                print(f"seg not spanning: {current_pos}, {seg}")
                return False
            current_pos = seg.end
        if chrom_nonforbidden_boundaries['end'] - current_pos >= allowance:
            print(f"last seg not spanning: {current_pos}, {chrom_nonforbidden_boundaries['end']}")
            return False
    return True

def batch_validate_OMKar_output(mk_dir, cont_allowance=50000, span_allowance=50000,
                                forbidden_region_file=get_metadata_file_path('acrocentric_telo_cen.bed')):
    def validate_single_file(filepath):
        path_list, segment_dict = read_OMKar_output(filepath, return_segment_dict=True)
        for index, segment in segment_dict.items():
            if not segment.direction():
                print(f"inverted direction: {segment}")
                return False
        if not check_continous(segment_dict, allowance=cont_allowance):
            return False
        if not check_spanning(segment_dict, forbidden_region_file, allowance=span_allowance):
            return False
        return True

    file_names = os.listdir(mk_dir)
    all_true = True
    false_file_paths = []
    for omkar_output in file_names:
        omkar_output_filepath = os.path.join(mk_dir, omkar_output)
        if not validate_single_file(omkar_output_filepath):
            all_true = False
            false_file_paths.append(omkar_output_filepath)
            print(f"^FALSE: {omkar_output}")
        else:
            print(f"^TRUE: {omkar_output}")

    print(f"***ALL_RETURNED_CORRECT: {all_true}")
    return false_file_paths

#################################POST-PROCESSING################################

def rotate_and_bin_path(path_list, forbidden_region_file=get_metadata_file_path('acrocentric_telo_cen.bed'), return_rotated_idx=False):
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

def post_process_OMKar_output(path_list, gap_allowance=5, isolate_centromere=True,
                              forbidden_region_file=get_metadata_file_path('acrocentric_telo_cen.bed')):
    ### invert paths that are output as backward
    tmp_path_list = []
    for path in path_list:
        tmp_path = path.duplicate()
        tmp_path_list.append(tmp_path)
    rotated_path_idx = rotate_and_bin_path(tmp_path_list, return_rotated_idx=True)
    for path_idx in rotated_path_idx:
        path = path_list[path_idx]
        path.reverse()

    ### merge segments that always appear together, in the same orientation
    ## enumerate longest contig/sublist of segments that are cont. 1) in orientation, 2) in positions (<= 5 bp gap) and Chr
    sublists = []
    for path in path_list:
        start_seg_idx = 0
        while start_seg_idx < len(path.linear_path.segments):
            end_seg_idx = legal_contig_extension(path.linear_path.segments, start_seg_idx, gap_allowance)
            sublist = path.linear_path.segments[start_seg_idx:end_seg_idx+1]
            ## reorient each sublist so that they are each in the + direction
            oriented_sublist = positively_orient_sublist(sublist)
            if oriented_sublist not in sublists:
                # prevent duplicates (this may not be an efficient implementation)
                sublists.append(oriented_sublist)
            start_seg_idx = end_seg_idx + 1
    sublists = sorted(sublists, key=lambda x: -len(x))
    ### break up all sublists if they pairwise intersect with another, until no further breaking is available
    ### then, remove all sublists with just one element
    while True:
        new_sublists_to_add = []
        sublist_to_remove_idx = set()
        ## all pairwise comparison
        for sublist_idx, sublist in enumerate(sublists):
            for compare_sublist_idx in range(sublist_idx+1, len(sublists)):
                compare_sublist = sublists[compare_sublist_idx]
                splitted_sublists = sublist_breaking(sublist, compare_sublist)
                if splitted_sublists:
                    sublist_to_remove_idx.add(sublist_idx)
                    sublist_to_remove_idx.add(compare_sublist_idx)
                    new_sublists_to_add += splitted_sublists

        sublist_to_remove_idx = list(sublist_to_remove_idx)
        sublist_to_remove_idx = sorted(sublist_to_remove_idx, reverse=True)
        if not sublist_to_remove_idx:
            break  # all intersections removed
        for sublist_to_remove_itr in sublist_to_remove_idx:
            sublists.pop(sublist_to_remove_itr)
        for sublist_itr in new_sublists_to_add:
            if sublist_itr not in sublists:
                sublists.append(sublist_itr)
        sublists = sorted(sublists, key=lambda x: -len(x))
    sublists_to_merge = [sub for sub in sublists if len(sub) >= 2]
    #
    # sublists_to_merge = []
    # while len(sublists) > 0:
    #     sublist = sublists.pop(0)
    #     if len(sublist) < 2:
    #         break  # all remainings are <2, so finding merge segments is completed
    #     sublist_to_be_removed = False
    #     sublist_to_compare_removal_idx = []
    #     sublist_to_add = []
    #     for compare_idx, sublist_to_compare in enumerate(sublists):
    #         # sublist_to_compare is <= in len() to sublist due to max-heap structure
    #
    #         if intersection_sublist:
    #             # remove the larger sublist
    #             sublist_to_be_removed = True
    #             if not is_substring:
    #                 # only keep the smaller sublist when it is a substring of the larger one, OW keep the intersection
    #                 sublist_to_compare_removal_idx.append(compare_idx)
    #                 sublist_to_add.append(intersection_sublist)
    #
    #     if not sublist_to_be_removed:
    #         sublists_to_merge.append(sublist)
    #     else:
    #         if sublist_to_compare_removal_idx:
    #             sublist_to_compare_removal_idx = sorted(sublist_to_compare_removal_idx, reverse=True)
    #             for idx in sublist_to_compare_removal_idx:
    #                 sublists.pop(idx)
    #             for sublist_itr in sublist_to_add:
    #                 sublists.append(sublist_itr)
    #             sublists = sorted(sublists, key=lambda x: -len(x))  # can be optimized

    ### find all occurences and merge
    # also scan for inverted segment
    inverted_sublists_to_merge = [invert_sublist(sublist) for sublist in sublists_to_merge]
    for path in path_list:
        ## find all ranges to be merged
        merging_ranges = []  # inclusive ranges
        for sublist in inverted_sublists_to_merge:
            inverted_starts = find_all_indices(path.linear_path.segments, sublist[0])
            for start in inverted_starts:
                # sanity check, disable for performance
                if path.linear_path.segments[start:start + len(sublist)] != sublist:
                    raise RuntimeError()
                merging_ranges.append((start, start + len(sublist) - 1))
        for sublist in sublists_to_merge:
            forward_starts = find_all_indices(path.linear_path.segments, sublist[0])
            for start in forward_starts:
                # sanity check, disable for performance
                if path.linear_path.segments[start:start + len(sublist)] != sublist:
                    segs_str = "".join([str(seg) for seg in path.linear_path.segments[start:start + len(sublist)]])
                    subl_str = "".join([str(seg) for seg in sublist])
                    print(f"segs: {segs_str}")
                    print(f"subl: {subl_str}")
                    raise RuntimeError()
                merging_ranges.append((start, start + len(sublist) - 1))

        ## merge
        merging_ranges = sorted(merging_ranges, key=lambda x: x[0], reverse=True)
        for merge_range in merging_ranges:
            # pop from the back to preserve list indices; ranges are also non-intersecting
            segments_to_merge = []
            for seg_idx in range(merge_range[1], merge_range[0]-1, -1):
                segments_to_merge.append(path.linear_path.segments.pop(seg_idx))
            segments_to_merge = segments_to_merge[::-1]
            merged_segment = merge_segments(segments_to_merge)
            path.linear_path.segments.insert(merge_range[0], merged_segment)

    ### create breakpoints at centromere boundaries
    if isolate_centromere:
        centromere_boundaries = get_centromere_boundaries(forbidden_region_file)
        for path in path_list:
            segment_arm = path.linear_path
            for chrom, cen_boundary in centromere_boundaries.items():
                segment_arm.introduce_breakpoint(chrom, cen_boundary['start'], 'start')
                segment_arm.introduce_breakpoint(chrom, cen_boundary['end'], 'end')

    ### collect segment_list
    all_segments = []
    for path in path_list:
        for seg in path.linear_path.segments:
            new_seg = seg.duplicate()
            if not new_seg.direction():
                new_seg.invert()
            if new_seg not in all_segments:
                all_segments.append(new_seg)
    all_segments = sorted(all_segments)
    segment_obj_to_idx_dict = {seg: idx + 1 for idx, seg in enumerate(all_segments)}
    return path_list, segment_obj_to_idx_dict

def batch_post_process_OMKar_output(omkar_output_dir, processing_output_dir, gap_merge_allowance=5,
                                    isolate_centromere=True,
                                    forbidden_region_file=get_metadata_file_path('acrocentric_telo_cen.bed')):
    os.makedirs(processing_output_dir, exist_ok=True)

    ### validate to skip files with issues
    files_with_issues = batch_validate_OMKar_output(omkar_output_dir)
    files_with_issues = [os.path.basename(file_name) for file_name in files_with_issues]

    for file_name in os.listdir(omkar_output_dir):
        if file_name in files_with_issues:
            print(f"skipping file with issue: {file_name}")
            continue
        print(f"post-processing: {file_name}")
        input_filepath = os.path.join(omkar_output_dir, file_name)
        output_filepath = os.path.join(processing_output_dir, file_name)
        path_list, segment_dict = read_OMKar_output(input_filepath, return_segment_dict=True)
        processed_path_list, segment_obj_to_idx_dict = post_process_OMKar_output(path_list,
                                                                                 gap_allowance=gap_merge_allowance,
                                                                                 isolate_centromere=isolate_centromere,
                                                                                 forbidden_region_file=forbidden_region_file)
        write_MK_file(output_filepath, processed_path_list, segment_obj_to_idx_dict)


def find_all_indices(lst, element):
    return [i for i, x in enumerate(lst) if x == element]

def legal_contig_extension(segments, start_idx, gap_allowance):
    c_seg = segments[start_idx]
    orientation = c_seg.direction()
    end_idx = start_idx
    for seg in segments[start_idx+1:]:
        if seg.direction() != orientation:
            break
        if orientation:
            if seg.start - c_seg.end > gap_allowance or seg.start - c_seg.end < 0:
                break
        else:
            if c_seg.end - seg.start > gap_allowance or c_seg.end - seg.start < 0:
                break
        c_seg = seg
        end_idx += 1
    return end_idx

def sublist_breaking(large_sublist, small_sublist):
    """
    @param large_sublist:
    @param small_sublist:
    @return: breakup the two sublists into multiple smaller sublists, so there is no intersection in the set(smaller sublists)
    If no intersection present from the start, return []
    """
    ## impossible to have intersection at two, non-contiguous locations, as both sublists are contigs
    intersection_start = -1
    intersection_end = -1
    for idx, seg in enumerate(large_sublist):
        if seg in small_sublist:
            intersection_start = idx
            break
    for idx, seg in enumerate(large_sublist[intersection_start:]):
        if seg not in small_sublist:
            intersection_end = intersection_start + idx  # not inclusive
            break
    if intersection_start != -1:
        if intersection_end == -1:
            intersection_sublist = large_sublist[intersection_start:]
        else:
            intersection_sublist = large_sublist[intersection_start:intersection_end]
    else:
        intersection_sublist = []
    if not intersection_sublist:
        return []

    small_sublist_start_idx = small_sublist.index(intersection_sublist[0])
    small_sublist_end_idx = small_sublist.index(intersection_sublist[-1]) + 1  # not inclusive

    splitted_sublists = [intersection_sublist]
    if small_sublist_start_idx != 0:
        new_sublist = small_sublist[:small_sublist_start_idx]
        if new_sublist not in splitted_sublists and len(new_sublist) >= 2:
            splitted_sublists.append(new_sublist)
    if small_sublist_end_idx != len(small_sublist):
        new_sublist = small_sublist[small_sublist_end_idx:]
        if new_sublist not in splitted_sublists and len(new_sublist) >= 2:
            splitted_sublists.append(new_sublist)
    if intersection_start != 0:
        new_sublist = large_sublist[:intersection_start]
        if new_sublist not in splitted_sublists and len(new_sublist) >= 2:
            splitted_sublists.append(new_sublist)
    if intersection_end != -1:
        new_sublist = large_sublist[intersection_end:]
        if new_sublist not in splitted_sublists and len(new_sublist) >= 2:
            splitted_sublists.append(new_sublist)

    return splitted_sublists

def invert_sublist(sublist):
    new_lst = []
    for seg in sublist[::-1]:
        new_seg = seg.duplicate()
        new_seg.invert()
        new_lst.append(new_seg)
    return new_lst

def positively_orient_sublist(sublist):
    ## assumes sublist is legal
    orientation = sublist[0].direction()
    if orientation:
        return sublist
    else:
        return invert_sublist(sublist)

def merge_segments(segment_sublist):
    def concatenate_kt_index(segments):
        return ''.join(segment.kt_index for segment in segments)

    ## assumes the gap/cont. were checked earlier
    return Segment(segment_sublist[0].chr_name,
                   segment_sublist[0].start,
                   segment_sublist[-1].end,
                   kt_index=concatenate_kt_index(segment_sublist))

####################################################################################

def write_MK_file(output_path, path_list, segment_obj_to_idx_dict):
    output_str = "Segment\tNumber\tChromosome\tStart\tEnd\tStartNode\tEndNode\n"
    for seg_obj, seg_idx in segment_obj_to_idx_dict.items():
        chrom = seg_obj.chr_name.replace("Chr", "")
        if chrom == "X":
            chrom = "23"
        elif chrom == "Y":
            chrom = "24"
        output_str += f"Segment\t{seg_idx}\t{chrom}\t{seg_obj.start}\t{seg_obj.end}\n"
    for idx, path in enumerate(path_list):
        seg_string = []
        for seg_obj in path.linear_path.segments:
            if not seg_obj.direction():
                new_seg = seg_obj.duplicate()
                new_seg.invert()
                seg_string.append(f"{segment_obj_to_idx_dict[new_seg]}-")
            else:
                seg_string.append(f"{segment_obj_to_idx_dict[seg_obj]}+")
        seg_string = " ".join(seg_string)
        output_str += f"Path{idx+1} = {seg_string}\n"
    output_str = output_str.strip()
    with open(output_path, "w") as fp_write:
        fp_write.write(output_str)

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


def read_OMKar_to_indexed_list(OMKar_output_file, forbidden_region_file=get_metadata_file_path('acrocentric_telo_cen.bed')):
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


################################TESTS############################


def test():
    path_list = read_OMKar_output("/media/zhaoyang-new/workspace/KarSim/KarComparator/new_data_files/OMKar_testbuild3/23X_15q26_overgrowth_r1.1.txt")
    rotate_and_bin_path(path_list, "Metadata/merged_forbidden_regions_unique.bed")
    # report_centromere_anomaly(path_list)
    for path in path_list:
        print(path)


def test_read_OMKar_output():
    path_list, segment_list = read_OMKar_output("sample_input/23Y_Cri_du_Chat_r1.1.txt", return_segment_dict=True)
    if batch_validate_OMKar_output(path_list, segment_dict=segment_list) == True:
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
    
