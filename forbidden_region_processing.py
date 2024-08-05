from .Structures import Arm, Segment, Path
from .utils import *


def read_forbidden_regions(forbidden_region_file) -> Arm:
    segment_list = []
    with open(forbidden_region_file) as fp_read:
        fp_read.readline()  # skip index line
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            new_segment = Segment(str(line[0]), int(line[1]), int(line[2]), str(line[3]))
            segment_list.append(new_segment)
    return Arm(segment_list, 'forbidden_regions')


def output_forbidden_regions_from_arm(input_arm: Arm, output_file_path):
    with open(output_file_path, 'w') as fp_write:
        fp_write.write('Chr\tStartPos\tEndPos\tType\n')
        for segment_itr in input_arm.segments:
            output_str = "{}\t{}\t{}\t{}\n".format(segment_itr.chr_name,
                                                   str(segment_itr.start),
                                                   str(segment_itr.end),
                                                   segment_itr.segment_type)
            fp_write.write(output_str)


def label_path_with_forbidden_regions(input_path_list: [Path], forbidden_region_file):
    forbidden_regions_path = Path(read_forbidden_regions(forbidden_region_file), 'forbidden_regions', 'forbidden_regions')
    for path_itr in input_path_list:
        # breakup into disjoint segments
        path_itr.generate_mutual_breakpoints(forbidden_regions_path, mutual=True)

        # label forbidden regions
        for path_segment_itr in path_itr.linear_path.segments:
            labeled = False
            for forbidden_region_segment_itr in forbidden_regions_path.linear_path.segments:
                if path_segment_itr.same_segment_ignore_dir(forbidden_region_segment_itr):
                    path_segment_itr.segment_type = forbidden_region_segment_itr.segment_type
                    labeled = True
            if not labeled:
                path_segment_itr.segment_type = 'arm_region'


def get_chr_length_from_forbidden_file(input_chr_name, forbidden_region_file=get_metadata_file_path("acrocentric_telo_cen.bed")):
    with open(forbidden_region_file) as fp_read:
        fp_read.readline()  # skip index line
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            chrom = line[0]
            end_pos = int(line[2])
            seg_type = line[3]
            if chrom.lower() == input_chr_name.lower():
                if seg_type == 'telomere2':
                    return end_pos
    return None


def get_prefix_suffix_forbidden_boundaries(forbidden_region_file=get_metadata_file_path("acrocentric_telo_cen.bed")):
    """
    used for creating source and sink nodes in comparison
    for prefix and suffix forbidden regions, start boundary is the last bp of the prefix forbidden segment; end boundary is the first bp of the suffix
    :param forbidden_region_file: assumes chr order starts with 1, and segment in increasing order; all the same chr segments are in the same block
    :return:
    """
    forbidden_region_arm = read_forbidden_regions(forbidden_region_file)
    boundaries = {f"Chr{i}": {'start': -1, 'end': -1} for i in list(range(1, 23)) + ['X', 'Y']}  # for every chromosome, there is a start + end boundary
    c_chr = 'Chr1'
    c_chr_start = False
    c_chr_end = False
    for seg in forbidden_region_arm.segments:
        seg_chr = seg.chr_name
        seg_type = seg.segment_type

        if seg_chr != c_chr:
            # chr changed
            c_chr = seg_chr
            c_chr_start = False
            c_chr_end = False

        # non-acrocentric start
        if (not c_chr_start) and seg_type.startswith('telomere'):
            boundaries[c_chr]['start'] = seg.end
            c_chr_start = True
            continue

        # acrocentric start
        if (not c_chr_start) and seg_type.startswith('acrocentric-telomere'):
            continue
        elif (not c_chr_start) and seg_type.startswith('acrocentric-centromere'):
            boundaries[c_chr]['start'] = seg.end
            c_chr_start = True
            continue

        # end
        if c_chr_start and (not c_chr_end) and seg_type.startswith('telomere'):
            boundaries[c_chr]['end'] = seg.start
            c_chr_end = True
            continue

    # check if all populated
    for boundary in boundaries.values():
        if boundary['start'] == -1 or boundary['end'] == -1:
            raise RuntimeError()
    return boundaries


def test():
    print(read_forbidden_regions(get_metadata_file_path("acrocentric_telo_cen.bed")))


if __name__ == "__main__":
    test()
