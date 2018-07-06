#!/usr/bin/env python

import os
import sys
import argparse
import re
import pysam
import misopy.sam_utils as sam_utils


def readRegionList(infile):
    exon_dict = {}
    with open(infile) as f:
        for line in f:
            F = line.rstrip().split("\t")
            region = [F[0]]
            region = region + re.split("[:-]", F[0])
            exon_dict[F[1]] = region
    return exon_dict


def getMate(samfile, read):
    if read.is_paired:
        mate = None
        pointer = samfile.tell()
        try:
            mate = samfile.mate(read)
        except ValueError:
            mate = None
        finally:
            samfile.seek(pointer)
        return mate
    else:
        return None


def isCD45RO(read, exon_dict):
    """ read: Exon3's reads jump to Exon7 """
    """ overlap() deprecated, use get_overlap() instead """
    nBaseE3 = read.get_overlap(int(exon_dict["PTPRC::Exon3"][2])-1, int(exon_dict["PTPRC::Exon3"][3]))
    nBaseE4 = read.get_overlap(int(exon_dict["PTPRC::Exon4"][2])-1, int(exon_dict["PTPRC::Exon4"][3]))
    nBaseE5 = read.get_overlap(int(exon_dict["PTPRC::Exon5"][2])-1, int(exon_dict["PTPRC::Exon5"][3]))
    nBaseE6 = read.get_overlap(int(exon_dict["PTPRC::Exon6"][2])-1, int(exon_dict["PTPRC::Exon6"][3]))
    nBaseE7 = read.get_overlap(int(exon_dict["PTPRC::Exon7"][2])-1, int(exon_dict["PTPRC::Exon7"][3]))
    if nBaseE3 > 0 and nBaseE4 == 0 and nBaseE5 == 0 and nBaseE6 == 0 and nBaseE7 > 0:
        return True
    else:
        return False


def isCD45RA(read, exon_dict):
    """ read: Exon3's reads jump to Exon4 """
    nBaseE3 = read.get_overlap(int(exon_dict["PTPRC::Exon3"][2])-1, int(exon_dict["PTPRC::Exon3"][3]))
    nBaseE4 = read.get_overlap(int(exon_dict["PTPRC::Exon4"][2])-1, int(exon_dict["PTPRC::Exon4"][3]))
    nBaseE5 = read.get_overlap(int(exon_dict["PTPRC::Exon5"][2])-1, int(exon_dict["PTPRC::Exon5"][3]))
    nBaseE6 = read.get_overlap(int(exon_dict["PTPRC::Exon6"][2])-1, int(exon_dict["PTPRC::Exon6"][3]))
    nBaseE7 = read.get_overlap(int(exon_dict["PTPRC::Exon7"][2])-1, int(exon_dict["PTPRC::Exon7"][3]))
    ####if nBaseE4 > 0 and nBaseE5 == 0 and nBaseE6 == 0 and nBaseE7 > 0:
    if nBaseE3 > 0 and nBaseE4 > 0:
        return True
    else:
        return False


def CD45StatusByJunction(infile, exon_dict, verbose, filter_reads):
    samfile = pysam.Samfile(infile, "rb")
    #it = samfile.fetch()
    it = samfile.fetch(region="chr1:198608098-198726605")
    n0 = 0
    nJunc = 0
    typeO = []
    typeA = []
    typeB_BC = []
    type56 = []
    type57 = []
    for read in it:
        if filter_reads:
            # Skip reads with no CIGAR string
            if read.cigar is None and verbose:
                print "#Skipping read with no CIGAR string: %s\t%s" %(read.cigar,read.qname)
                continue
            cigar_str = sam_utils.sam_cigar_to_str(read.cigar)
            if ("N" in cigar_str) and (cigar_str.count("N") > 1) and verbose:
                print "#Skipping read with multiple junctions crossed: %s\t%s" %(cigar_str,read.qname)
                continue
            # Check if the read contains an insertion (I)
            # or deletion (D) -- if so, skip it
            for cigar_part in read.cigar:
                if (cigar_part[0] == 1 or cigar_part[1] == 2) and verbose:
                    print "#Skipping read with CIGAR %s\t%s" %(cigar_str,read.qname)
                    continue

        n0 += 1
        ##isTypeO = isCD45RO(read, exon_dict)
        ##isTypeA = isCD45RA(read, exon_dict)
        ##typeO.append(isTypeO)
        ##typeA.append(isTypeA)
        nBaseE3 = read.get_overlap(int(exon_dict["PTPRC::Exon3"][2])-1, int(exon_dict["PTPRC::Exon3"][3]))
        nBaseE4 = read.get_overlap(int(exon_dict["PTPRC::Exon4"][2])-1, int(exon_dict["PTPRC::Exon4"][3]))
        nBaseE5 = read.get_overlap(int(exon_dict["PTPRC::Exon5"][2])-1, int(exon_dict["PTPRC::Exon5"][3]))
        nBaseE6 = read.get_overlap(int(exon_dict["PTPRC::Exon6"][2])-1, int(exon_dict["PTPRC::Exon6"][3]))
        nBaseE7 = read.get_overlap(int(exon_dict["PTPRC::Exon7"][2])-1, int(exon_dict["PTPRC::Exon7"][3]))
        if sum([i>0 for i in (nBaseE3,nBaseE4,nBaseE5,nBaseE6,nBaseE7)]) > 1:
            nJunc += 1
        haveRO  = False
        haveRA  = False
        haveRB_BC = False
        have57 = False
        ##haveRBC = False
        ##haveRB  = False
        """ RO: Exon3's reads jump to Exon7 """
        if nBaseE3 > 0 and nBaseE4 == 0 and nBaseE5 == 0 and nBaseE6 == 0 and nBaseE7 > 0:
            haveRO = True
        else:
            haveRO = False
        """ RA: Exon3's reads jump to Exon4 """
        if nBaseE3 > 0 and nBaseE4 > 0:
            haveRA = True
        else:
            haveRA = False
        """ RB_BC: Exon3's reads jump to Exon5 """
        if nBaseE3 > 0 and nBaseE4 == 0 and nBaseE5 > 0:
            haveRB_BC = True
        else:
            haveRB_BC = False
        """ 56: Exon5's reads jump to Exon6 """
        if nBaseE5 > 0 and nBaseE6 > 0:
            have56 = True
        else:
            have56 = False
        """ 57: Exon5's reads jump to Exon7 """
        if nBaseE5 > 0 and nBaseE6 == 0 and nBaseE7 > 0:
            have57 = True
        else:
            have57 = False
        typeO.append(haveRO)
        typeA.append(haveRA)
        typeB_BC.append(haveRB_BC)
        type56.append(have56)
        type57.append(have57)
        if verbose:
            print "#%s\t%s\t%s\t%s\t%s\t%s" % (read.qname, str(haveRO), str(haveRA), str(haveRB_BC),str(have56),str(have57))
            #print >> sys.stderr, mread
    samfile.close()
    if True in typeO:
        return ("CD45RO", sum(typeO), sum(typeA), sum(typeB_BC), sum(type56), sum(type57), n0, nJunc)
    elif True in typeA:
        return ("CD45RA", sum(typeO), sum(typeA), sum(typeB_BC), sum(type56), sum(type57), n0, nJunc)
    elif True in typeB_BC:
        if True in type57 and True not in type56:
            return ("CD45RB", sum(typeO), sum(typeA), sum(typeB_BC), sum(type56), sum(type57), n0, nJunc)
        elif True not in type57 and True in type56:
            return ("CD45RBC", sum(typeO), sum(typeA), sum(typeB_BC), sum(type56), sum(type57), n0, nJunc)
    ## not enough junction reads ?
    if nJunc == 0:
        return ("lowCov", sum(typeO), sum(typeA), sum(typeB_BC), sum(type56), sum(type57), n0, nJunc)
    else:
        return ("unknown", sum(typeO), sum(typeA), sum(typeB_BC), sum(type56), sum(type57), n0, nJunc)

def main():
    parser = argparse.ArgumentParser(description='Type the CD45 status.')
    parser.add_argument("--sample", dest="sample_id", type=str, nargs='?',
                        default="SAMPLE", help="sample id")
    parser.add_argument("--input-bam", dest="input_bam", type=str, nargs='?',
                        required=True, default=None, help="Input bam file.")
    parser.add_argument("--region-list", dest="region_list", type=str, nargs='?',
                        required=True, default=None,
                        help="Region list, each line:\n\"chr1:198661476-198661502    PTPRC::Exon3\"")
    parser.add_argument("--filterReads", dest="filter_reads", action='store_true', default=False,
                        help="filter reads")
    parser.add_argument("--verbose", dest="verbose", action='store_true', default=False,
                        help="Verbose output")
    args = parser.parse_args()
    exon_dict = readRegionList(args.region_list)
    otype = CD45StatusByJunction(args.input_bam, exon_dict, args.verbose, args.filter_reads)
    print "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (args.sample_id, otype[0], otype[1], otype[2], otype[3], otype[4], otype[5], otype[6], otype[7])


if __name__ == '__main__':
    main()
