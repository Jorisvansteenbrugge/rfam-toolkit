#!/usr/bin/python
'''
Created on 7 Dec 2015

Description: A tool designed to support Rfam scanning both locally 
             as well as utilizing IBM's LSF platform.

@author: ikalvari
'''

import os
import sys
import subprocess

#-------------------------------------------------------------------------


LSF = False  # True: LSF execution, False: Local execution
CPU = 4  # CPUs to reserve on exec host
MEM = 4000  # Memory to be reserved
TMP_MEM = 9000  # TMP memory to reserve on LSF - 9GB

TMP_PATH = "/tmp"  # tmp directory path
POST_EXEC = "rm %s/$LSB_JOBID.*" % TMP_PATH  # rm /tmp/jobid.*"


CMSEARCH_EXEC = None
CMSCAN_EXEC = None

# Setting execution paths
if (LSF is True):
    # LSF exec path
    CMSEARCH_EXEC = "/path/to/lsf/cmsearch"
    CMSCAN_EXEC = "/path/to/lsf/cmscan"
else:
    # local Infernal exec path
    CMSEARCH_EXEC = "/path/to/infernal-1.1/src/cmsearch"
    CMSCAN_EXEC = "/path/to/infernal-1.1/src/cmscan"

CMSEARCH = 'CS'
CMSCAN = 'CC'


CMD_TEMPLATE = "%s --tblout %s.out --notextw --cut_ga --nohmmonly --rfam --cpu %d %s %s"

# LSF SPECIFIC
LSB_JOBID = "$LSB_JOBID"
GROUP_NAME = "/group_name"  # LSF Group Name, create in advance (e.g. /rfam_rnac)


# TO DO
PARAMS = None
BSUB_PARAMS = None

#----------------------------------FUNCTIONS------------------------------


def batch_search(cm_dir, seq_dir, multi=False, method=CMSEARCH, out_dir=None):
    '''
        This function runs the Infernal's CMSEARCH and CMSCAN commands on a 
        set of CMs (.cm) provided in cm_dir and a set of sequence files provided 
        in seq_dir (.fa) and places the output in the destination provided as 
        out_dir. params are not currently used..

        By default cmsearch is the search method, unless set otherwise
        cm_dir:     The path to the directory containing the cm model files
        seq_dir:    The path to the directory containing the fasta files 
        multi:      If set to True, the code works on multiple seq directories
        method:     Infernal's search method CMSEARCH or CMSCAN
        out_dir:    This should be the path of the output directory. If not 
                    provided, seq_dir is used instead
    '''

    # List available cms
    cms = os.listdir(cm_dir)

    # Filtering out cmscan binaries
    # if(method==CMSCAN): #perform as a check or only for CMSCAN ?
    cms = filter(lambda x: x.endswith(".cm"), cms)

    seqs = os.listdir(seq_dir)

    # script_dir_path = None
    tool = None

    # setting the tool
    if(method == CMSEARCH):
        tool = CMSEARCH_EXEC
    else:
        tool = CMSCAN_EXEC
        # perhaps call cm file filtering here

    for cm in cms:
        # get the path for a specific model file
        cm_path = os.path.join(cm_dir, cm)

        # family accession used for naming purposes
        fam_acc = cm.partition('.')[0]
        fam_outdir = None

        # create family output directory
        if(out_dir != None):
            fam_outdir = os.path.join(out_dir, fam_acc)
        else:
            fam_outdir = os.path.join(seq_dir, fam_acc)

        if(os.path.exists(fam_outdir) is False):
            try:
                os.mkdir(fam_outdir)
            except:
                # TO DO
                pass

        # in the case of multi seq_file will be a directory
        for seq_entity in seqs:

            if(multi is True):
                proj_path = os.path.join(seq_dir, seq_entity)
                seq_files = os.listdir(proj_path)
                proj_dest_path = os.path.join(fam_outdir, seq_entity)

                # create out_dir
                if(os.path.exists(proj_dest_path) is False):
                    try:
                        os.mkdir(proj_dest_path)
                    except:
                        # TO DO
                        pass

                for seq_file in seq_files:

                    # get the path for a specific sequence file
                    seq_file_path = os.path.join(proj_path, seq_file)
                    seq_filename = seq_file.partition('.')[0]

                    cmd = None

                    # run rfam_scanner on the LSF platform
                    if (LSF is True):
                        cmd = CMD_TEMPLATE % (tool, os.path.join(TMP_PATH, LSB_JOBID), CPU,
                                              os.path.join(
                                                  TMP_PATH, LSB_JOBID + ".cm"),
                                              os.path.join(TMP_PATH, LSB_JOBID + ".fa"))
                        sh_path = generate_job_script(
                            fam_acc + '_' + 
                            seq_filename, cm_path, seq_file_path,
                            cmd, out_dir=proj_dest_path)

                        cmd = None
                        cmd = "bsub < %s" % (sh_path)

                    # run rfam_scanner locally
                    else:
                        cmd = CMD_TEMPLATE % (tool, os.path.join(
                            proj_dest_path, fam_acc + '_' + seq_filename),
                            CPU, cm_path, seq_file_path)

                    subprocess.call(cmd, shell=True)  # local installation

                    cmd = None

            else:

                # get the path for a specific sequence file
                seq_file_path = os.path.join(seq_dir, seq_entity)
                seq_filename = seq_entity.partition('.')[0]

                # Call cmd generator here - TO DO

                cmd = None

                # run rfam_scanner on the LSF platform
                if (LSF is True):
                    cmd = CMD_TEMPLATE % (tool, os.path.join(TMP_PATH, LSB_JOBID), CPU,
                                          os.path.join(
                                              TMP_PATH, LSB_JOBID + ".cm"),
                                          os.path.join(TMP_PATH, LSB_JOBID + ".fa"))
                    sh_path = generate_job_script(
                        fam_acc + '_' + seq_filename, cm_path, seq_file_path, cmd, out_dir=fam_outdir)

                    cmd = None
                    cmd = "bsub < %s" % (sh_path)

                # run rfam_scanner locally
                else:
                    cmd = CMD_TEMPLATE % (tool, os.path.join(
                        fam_outdir, fam_acc + '_' + seq_filename), CPU, cm_path, seq_file_path)

                subprocess.call(cmd, shell=True)  # local installation

                cmd = None

#-------------------------------------------------------------------------


def generate_job_script(filename, cm_file_path, seq_file_path, rfam_cmd, out_dir=None):

    # Modify this to accept the bsub parameters
    '''
        This function generates the .sh script for a job to be submitted to the LSF Platform.

        filename:      A string specifying the file's name
        cm_file_path:  The path to a specific CM file
        seq_file_path: The path to a fasta sequence file
        rfam_cmd:      A string of the predefined Infernal command, modified to run Infernal with
                       the input
        out_dir:       The path to the output directory. If not provided, sequence file directory
                       is used instead
    '''

    filepath = os.path.join(out_dir, filename + ".sh")

    fp = open(filepath, 'w')

    fp.write("#!/bin/csh\n")
    fp.write("#BSUB -M %d\n" % (MEM))
    fp.write("#BSUB -R \"rusage[mem=%d,tmp=%d]\"\n" % (MEM, TMP_MEM))
    # copying .cm file to exec host
    fp.write("#BSUB -f \"%s > %s/%sJ.cm\"\n" % 
             (cm_file_path, TMP_PATH, chr(37)))
    # copying .fa file to exec host
    fp.write("#BSUB -f \"%s > %s/%sJ.fa\"\n" % 
             (seq_file_path, TMP_PATH, chr(37)))
    fp.write("#BSUB -f \"%s/%s.out < %s/%sJ.out\"\n" % 
             (out_dir, filename, TMP_PATH, chr(37)))  # copying .out file to outdir
    fp.write("#BSUB -f \"%s/%s.err < %s/%sJ.err\"\n" % 
             (out_dir, filename, TMP_PATH, chr(37)))  # copy .err file to outdir
    fp.write("#BSUB -e \"%s/%sJ.err\"\n" % (TMP_PATH, chr(37)))
    fp.write("#BSUB -Ep \"%s\"\n" % (POST_EXEC))  # erase inputs from /tmp
    fp.write("#BSUB -n %d\n" % (CPU))
    # moving this to bsub command
    fp.write("#BSUB -g \"%s\"\n\n" % (GROUP_NAME))
    fp.write(rfam_cmd)  # construct the command here

    fp.close()

    return filepath

#-------------------------------------------------------------------------


def build_infernal_cmd(cm_file, seq_file, outdir, method=CMSEARCH, params=PARAMS):
    '''
        TO DO
    '''

    search_cmd = ""

    # IMPLEMENT HERE

    return search_cmd

#-------------------------------------------------------------------------


def rfam_scanner_help():
    '''
        TO DO
    '''

    pass

#-------------------------------------------------------------------------


def main(cm_dir, proj_dir, multi, method, out_dir):  # Remove?

    batch_search(cm_dir, proj_dir, multi=multi, method=method, out_dir=out_dir)

#-------------------------------------------------------------------------

if __name__ == '__main__':

    multi = False
    method = CMSEARCH

    cm_dir = sys.argv[1]
    proj_dir = sys.argv[2]
    out_dir = sys.argv[3]

    if (len(sys.argv) == 5 or len(sys.argv) == 6):
        if (sys.argv[4] == '-m'):
            multi = True
    elif (len(sys.argv) == 2 and sys.argv[1] == "-h"):
        print "Help currently unavailable"
        # rfam_scanner_help()

    batch_search(cm_dir, proj_dir, multi=multi, method=method, out_dir=out_dir)
