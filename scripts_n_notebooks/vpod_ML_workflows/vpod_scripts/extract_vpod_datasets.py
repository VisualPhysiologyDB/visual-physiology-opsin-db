"""
extract_vpod_datasets.py : routines for formatting vpod datasets for training ML models

extract_vpod_datasets(mydb,
                           ignore_filter=False,
                           only_visual_opsins=True,
                           keep_conflicts=False,
                           output_dir_prefix='vpod_1.2_data_splits',
                           bovine_seq_str=">Bovine\nMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA\n",
                           invert_seq_str=">Squid\nMGRDLRDNETWWYNPSIVVHPHWREFDQVPDAVYYSLGIFIGICGIIGCGGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFLKKWIFGFAACKVYGFIGGIFGFMSIMTMAMISIDRYNVIGRPMAASKKMSHRRAFIMIIFVWLWSVLWAIGPIFGWGAYTLEGVLCNCSFDYISRDSTTRSNILCMFILGFFGPILIIFFCYFNIVMSVSNHEKEMAAMAKRLNAKELRKAQAGANAEMRLAKISIVIVSQFLLSWSPYAVVALLAQFGPLEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFREAISQTFPWVLTCCQFDDKETEDDKDAETEIPAGESSDAAPSADAAQMKEMMAMMQKMQQQQAAYPPQGYAPPPQGYPPQGYPPQGYPPQGYPPQGYPPPPQGAPPQGAPPAAPPQGVDNQAYQA\n",
                           meta_first_line_str="Seq_Id\tLambda_Max\tSpecies\tOpsin_Family\tPhylum\tClass\tAccession\tMutations\tProtein\tRefId\nBovine\t500.0000\tBos_tarus\tRh1\tChordata\tMammalia\tNM_001014890\t\tMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA\n",
                           invert_first_line_str="Seq_Id\tLambda_Max\tSpecies\tOpsin_Family\tPhylum\tClass\tAccession\tMutations\tProtein\tRefId\nSquid\t473.0000\tTodarodes_pacificus\tRh1\tMollusca\tCephalopoda\t\tX70498\tMGRDLRDNETWWYNPSIVVHPHWREFDQVPDAVYYSLGIFIGICGIIGCGGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFLKKWIFGFAACKVYGFIGGIFGFMSIMTMAMISIDRYNVIGRPMAASKKMSHRRAFIMIIFVWLWSVLWAIGPIFGWGAYTLEGVLCNCSFDYISRDSTTRSNILCMFILGFFGPILIIFFCYFNIVMSVSNHEKEMAAMAKRLNAKELRKAQAGANAEMRLAKISIVIVSQFLLSWSPYAVVALLAQFGPLEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFREAISQTFPWVLTCCQFDDKETEDDKDAETEIPAGESSDAAPSADAAQMKEMMAMMQKMQQQQAAYPPQGYAPPPQGYPPQGYPPQGYPPQGYPPQGYPPPPQGAPPQGAPPAAPPQGVDNQAYQA\n")
                           - Extracts and splits VPOD datasets based on various criteria and filters.
"""

import datetime
import os
import re
import mysql.connector
import fileinput
import shutil

def extract_vpod_datasets(mydb,
                           ignore_filter=False,
                           only_visual_opsins=True,
                           keep_conflicts=False,
                           output_dir_prefix='vpod_1.2_data_splits',
                           bovine_seq_str=">Bovine\nMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA\n",
                           invert_seq_str=">Squid\nMGRDLRDNETWWYNPSIVVHPHWREFDQVPDAVYYSLGIFIGICGIIGCGGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFLKKWIFGFAACKVYGFIGGIFGFMSIMTMAMISIDRYNVIGRPMAASKKMSHRRAFIMIIFVWLWSVLWAIGPIFGWGAYTLEGVLCNCSFDYISRDSTTRSNILCMFILGFFGPILIIFFCYFNIVMSVSNHEKEMAAMAKRLNAKELRKAQAGANAEMRLAKISIVIVSQFLLSWSPYAVVALLAQFGPLEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFREAISQTFPWVLTCCQFDDKETEDDKDAETEIPAGESSDAAPSADAAQMKEMMAMMQKMQQQQAAYPPQGYAPPPQGYPPQGYPPQGYPPQGYPPQGYPPPPQGAPPQGAPPAAPPQGVDNQAYQA\n",
                           meta_first_line_str="Seq_Id\tLambda_Max\tSpecies\tOpsin_Family\tPhylum\tClass\tAccession\tMutations\tProtein\tRefId\nBovine\t500.0000\tBos_tarus\tRh1\tChordata\tMammalia\tNM_001014890\t\tMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA\n",
                           invert_first_line_str="Seq_Id\tLambda_Max\tSpecies\tOpsin_Family\tPhylum\tClass\tAccession\tMutations\tProtein\tRefId\nSquid\t473.0000\tTodarodes_pacificus\tRh1\tMollusca\tCephalopoda\t\tX70498\tMGRDLRDNETWWYNPSIVVHPHWREFDQVPDAVYYSLGIFIGICGIIGCGGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFLKKWIFGFAACKVYGFIGGIFGFMSIMTMAMISIDRYNVIGRPMAASKKMSHRRAFIMIIFVWLWSVLWAIGPIFGWGAYTLEGVLCNCSFDYISRDSTTRSNILCMFILGFFGPILIIFFCYFNIVMSVSNHEKEMAAMAKRLNAKELRKAQAGANAEMRLAKISIVIVSQFLLSWSPYAVVALLAQFGPLEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFREAISQTFPWVLTCCQFDDKETEDDKDAETEIPAGESSDAAPSADAAQMKEMMAMMQKMQQQQAAYPPQGYAPPPQGYPPQGYPPQGYPPQGYPPQGYPPPPQGAPPQGAPPAAPPQGVDNQAYQA\n"):
    """
    Extracts and splits VPOD datasets based on various criteria and filters.

    Args:
        mydb: MySQL database connection object.
        ignore_filter (bool, optional): If True, redundant opsin data (identical sequences) will be kept. Defaults to False.
        only_visual_opsins (bool, optional): If True, only visual opsins will be kept, and non-visual opsins will be filtered out. Defaults to True.
        keep_conflicts (bool, optional): If True, redundant sequences with conflicting lmax values will be kept. Defaults to False.
        output_dir_prefix (str, optional): Prefix for the output directory name. Defaults to 'vpod_1.2_data_splits'.
        bovine_seq_str (str, optional): Bovine sequence string to write to output files. Defaults to default bovine sequence.
        invert_seq_str (str, optional): Invert sequence string to write to output files. Defaults to default squid sequence.
        meta_first_line_str (str, optional): Header string for most metadata files. Defaults to default meta header.
        invert_first_line_str (str, optional): Header string for invertebrate metadata files. Defaults to default invert header.

    Returns:
        str: Path to the directory where the output files are saved.
    """

    # directory preperation
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    seq_report_dir = str(f'{output_dir_prefix}_{dt_label}')
    os.makedirs(seq_report_dir, exist_ok=True)  # Added exist_ok=True to prevent error if directory exists

    # declaring all variables for different sequence data subsets and their metadata
    wd_output = f'{seq_report_dir}/wds.txt'
    sws_output = f'{seq_report_dir}/uss.txt'
    mws_output = f'{seq_report_dir}/mls.txt'
    rod_output = f'{seq_report_dir}/rod.txt'
    wd_ni_output = f'{seq_report_dir}/vert.txt'
    wt_vert_output = f'{seq_report_dir}/wt_vert.txt'
    inv_output = f'{seq_report_dir}/inv.txt'
    wt_inv_output = f'{seq_report_dir}/wt_inv.txt' # Not used in data_split_list, might be redundant
    nmoc_output = f'{seq_report_dir}/wt.txt'
    mut_output = f'{seq_report_dir}/mut_only.txt'
    wh_metadata = f'{seq_report_dir}/wds_meta.tsv'
    sw_metadata = f'{seq_report_dir}/uss_meta.tsv'
    mw_metadata = f'{seq_report_dir}/mls_meta.tsv'
    rh_metadata = f'{seq_report_dir}/rod_meta.tsv'
    wd_ni_metadata = f'{seq_report_dir}/vert_meta.tsv'
    inv_metadata = f'{seq_report_dir}/inv_meta.tsv'
    nmoc_metadata = f'{seq_report_dir}/wt_meta.tsv'
    mut_metadata = f'{seq_report_dir}/mut_meta.tsv'
    wt_vert_metadata = f'{seq_report_dir}/wt_vert_meta.tsv'
    redundant_datapoints = f'{seq_report_dir}/redundant_datapoints_log.tsv'


    data_split_list = [wd_output, sws_output, mws_output, rod_output, wd_ni_output, inv_output, nmoc_output, wt_vert_output, mut_output]
    meta_data_list = [wh_metadata, sw_metadata, mw_metadata, rh_metadata, wd_ni_metadata, inv_metadata, nmoc_metadata, wt_vert_metadata, mut_metadata]
    mnm_meta_list = [wh_metadata, nmoc_metadata, wd_ni_metadata, wt_vert_metadata, inv_metadata] # wt_vert_output should likely be wt_vert_metadata in meta lists
    mnm_meta_shorthand = ['wds', 'wt', 'vert', 'wt_vert', 'inv']
    # Setting the names for the headers at the top of each metadata file
    meta_first_line = meta_first_line_str
    invert_first_line = invert_first_line_str
    bovine_seq = bovine_seq_str

    # misc int variales for the counting loops
    m = 0
    s = 0
    l = 0
    r = 0
    c = 0
    z = 0
    q = 0
    wt_vert = 0
    mut = 0

    # regular expressions for filtering different opsin family types
    rod = re.compile('^Rh$|Rh[0-2]|exoRh')
    d = re.compile("^NM_001014890.2$|^NM_001014890$|^P02699.1$")
    sws_reg = re.compile('^SWS')
    uvs_reg = re.compile('^UVS')
    mws_reg = re.compile('^MWS')
    lws_reg = re.compile('^LWS')
    rh1_reg = re.compile('^Rh$|^Rh[0-1]|exoRh')
    rh2_reg = re.compile('^Rh2')
    non_viz = re.compile('^POps|MOps|GoC|PeroOps')

    sql = """
    SELECT DISTINCT
        o.genus, o.species, o.genefamily, o.accession,
        h.lamdamax, o.aa, o.phylum, o.class, h.mutations, o.refid
    FROM
        opsins o
    JOIN
        heterologous h
    ON
        o.accession = h.accession AND o.refid = h.refid;
    """
    mycursor = mydb.cursor()
    mycursor.execute(sql)
    myresult = mycursor.fetchall()
    seq_list = []
    lmax_list = []
    species_list = []
    acc_list = []
    gene_family_list = []
    conflicting_lmax = {}

    for x in myresult:
        skip = False
        lmax = x[4]
        try:
            index = seq_list.index(x[5].strip().replace(' ', ''))
            lmax_of_index = lmax_list[index]
            species_of_index = species_list[index]
            acc_of_index = acc_list[index]
            gene_family_of_index = gene_family_list[index]
        except ValueError: # Changed from except: to except ValueError: for clarity
            lmax_of_index = None

        species_check = str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '')
        if (x[5].strip().replace(' ', '') in seq_list) and (float(x[4]) != float(lmax_of_index)) and (lmax_of_index != None):
            # Redundant sequence found
            # Handle conflicting lmax values
            skip = True
            with open(redundant_datapoints, 'a+') as f:
                f.write(f'''Accessions {acc_of_index} and {x[3]} encode the same sequence but have conflicting lmax values ({lmax_of_index} and {lmax} respectively).\n''')
                if species_of_index != species_check:
                    f.write(f'''The two sequences are from different species ({species_of_index} and {species_check}).\n''')
                else:
                    f.write(f'''The two sequences are from the same species ({species_of_index}).\n''')
                f.write(f'''We will be taking the average of all conflicting values, and dropping this entry from the final datasplits.\n''')
                f.write(f'''Note - the gene families of the two sequences are {gene_family_of_index} and {x[2]} - if this raises concern, please double check the entries to resolve this issue.\n\n''')

            if not keep_conflicts: # Changed to if not keep_conflicts for readability
                # Store conflicting lmax values to modify original with average of all conflicting values
                if acc_of_index not in conflicting_lmax:
                    conflicting_lmax[acc_of_index] = [float(lmax_of_index), float(x[4])]  # Convert x[4] to float
                else:
                    conflicting_lmax[acc_of_index].append(float(x[4]))
            else:
                skip = False

        elif (x[5].strip().replace(' ', '') in seq_list) and (x[4] == lmax_of_index):
            # Redundant sequence found
            # Handle identical sequences and lmax values
            # Skipping exactly redundant entries
            skip = True
            with open(redundant_datapoints, 'a+') as f:
                f.write(f'''Accessions {acc_of_index} and {x[3]} encode the same sequence and have the same lmax values ({lmax_of_index}).''')
                if species_of_index != species_check:
                    f.write(f'''\nDespite that the two sequences are from different species ({species_of_index} and {species_check}) we will skip this entry but take the value from the first entry.\nNote - the gene families of the two sequences are {gene_family_of_index} and {x[2]} -  if this raises concern, please double check the entries to resolve this issue.\n\n''')
                else:
                    f.write(f'''Because the two sequences are from the same species ({species_of_index}) we will skip this entry and only keep the value from the first entry.\n\n''')
        else:
            # Store new sequence data
            species_list.append(species_check)
            gene_family_list.append(x[2])
            acc_list.append(x[3])
            lmax_list.append(lmax)
            seq_list.append(x[5].strip().replace(' ', ''))

        if ignore_filter: # Changed to if ignore_filter for readability
            skip = False

        if only_visual_opsins: # Changed to if only_visual_opsins for readability
            if non_viz.match(x[2]):
                skip = True

        if (lmax == 0) or skip: # Changed to if (lmax == 0) or skip for readability
            pass
        else:
            # SEQUENCE-DATA SECTION
            # whole-dataset
            with open(wd_output, 'a') as f:
                if m == 0:
                    f.write(bovine_seq)
                if (d.match(x[3])):
                    pass
                else:
                    m += 1
                    # This makes the fasta format file
                    seq = ">S" + str(m)
                    f.write(seq)
                    seq2 = str('\n' + x[5] + '\n')
                    f.write(seq2)

            # vertebrate dataset
            with open(wd_ni_output, 'a') as f:
                if x[6] != "Chordata" or d.match(x[3]):
                    pass
                else:
                    if c == 0:
                        f.write(bovine_seq)
                    c += 1
                    # This makes the fasta format file
                    seq = ">S" + str(c)
                    f.write(seq)
                    seq2 = str('\n' + x[5] + '\n')
                    f.write(seq2)

            # invertebrate dataset
            with open(inv_output, 'a') as f:
                if (x[6] != "Chordata"):
                    if q == 0:
                        f.write(invert_seq_str)
                    q += 1
                    # This makes the fasta format file
                    seq = ">S" + str(q)
                    f.write(seq)
                    seq2 = str('\n' + x[5] + '\n')
                    f.write(seq2)
                else:
                    pass

            # wild-type dataset
            with open(nmoc_output, 'a') as f:
                if ((len(x[8]) > 1) or ('-' in x[3] and '_(' not in x[3])):
                    pass
                else:
                    if z == 0:
                        f.write(bovine_seq)
                    if (d.match(x[3])):
                        pass
                    else:
                        z += 1
                        # This makes the fasta format file
                        seq = ">S" + str(z)
                        f.write(seq)
                        seq2 = str('\n' + x[5] + '\n')
                        f.write(seq2)

            # wild-type vertebrates dataset
            with open(wt_vert_output, 'a') as f:
                if (x[6] != "Chordata") or d.match(x[3]) or ((len(x[8]) > 1) or ('-' in x[3] and '_(' not in x[3])):
                    pass
                else:
                    if wt_vert == 0:
                        f.write(bovine_seq)
                    if (d.match(x[3])):
                        pass
                    else:
                        wt_vert += 1
                        # This makes the fasta format file
                        seq = ">S" + str(wt_vert)
                        f.write(seq)
                        seq2 = str('\n' + x[5] + '\n')
                        f.write(seq2)

            # just mutants dataset
            with open(mut_output, 'a') as f:
                if ((len(x[8]) > 1) or ('-' in x[3] and '_(' not in x[3])):
                    if mut == 0:
                        f.write(bovine_seq)
                    if (d.match(x[3])):
                        pass
                    else:
                        mut += 1
                        # This makes the fasta format file
                        seq = ">M" + str(mut)
                        f.write(seq)
                        seq2 = str('\n' + x[5] + '\n')
                        f.write(seq2)
                else:
                    pass

            # UVS and SWS dataset
            with open(sws_output, 'a') as f:
                p = re.compile('^SWS|^UVS')
                if p.match(x[2]) and x[6] == "Chordata":
                    s += 1
                    if s == 1:
                        f.write(bovine_seq)
                    # This makes the fasta format file
                    seq = ">S" + str(s)
                    f.write(seq)
                    seq2 = str('\n' + x[5] + '\n')
                    f.write(seq2)

            # MWS and LWS dataset
            with open(mws_output, 'a') as f:
                p = re.compile('^MWS|^LWS')
                if p.match(x[2]) and x[6] == "Chordata":
                    l += 1
                    if l == 1:
                        f.write(bovine_seq)
                    # This makes the fasta format file
                    seq = ">S" + str(l)
                    f.write(seq)
                    seq2 = str('\n' + x[5] + '\n')
                    f.write(seq2)

            # rods dataset
            with open(rod_output, 'a') as f:
                p = re.compile('Rh[0-2]|exoRh')
                if p.match(x[2]):
                    if r == 0:
                        f.write(bovine_seq)
                    if (x[6] != "Chordata") or d.match(x[3]):
                        pass
                    else:
                        r += 1
                        # This makes the fasta format file
                        seq = ">S" + str(r)
                        f.write(seq)
                        seq2 = str('\n' + x[5] + '\n')
                        f.write(seq2)

            # METADATA SECTION - same idea and naming convention for the files as above.
            with open(wh_metadata, 'a') as g:
                if m == 1:
                    g.write(meta_first_line)
                if (d.match(x[3])):
                    pass
                else:
                    md = str("S" + str(m) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                    g.write(md)

            with open(wd_ni_metadata, 'a') as g:
                if x[6] != "Chordata" or d.match(x[3]):
                    pass
                else:
                    if c == 1:
                        g.write(meta_first_line)

                    md = str("S" + str(c) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                    g.write(md)

            with open(inv_metadata, 'a') as g:
                if (x[6] != "Chordata"):
                    if q == 1:
                        g.write(invert_first_line)
                    md = str("S" + str(q) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                    g.write(md)
                else:
                    pass

            with open(sw_metadata, 'a') as g:
                p = re.compile('^SWS|^UVS')
                if p.match(x[2]) and x[6] == "Chordata":
                    if s == 1:
                        g.write(meta_first_line)
                    md = str("S" + str(s) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                    g.write(md)

            with open(mw_metadata, 'a') as g:
                p = re.compile('^MWS|^LWS')
                if p.match(x[2]) and x[6] == "Chordata":
                    if l == 1:
                        g.write(meta_first_line)
                    md = str("S" + str(l) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                    g.write(md)

            with open(rh_metadata, 'a') as g:
                p = re.compile('Rh[0-3]|exoRh')

                if p.match(x[2]):
                    if r == 1:
                        g.write(meta_first_line)
                    if (x[6] != "Chordata") or d.match(x[3]):
                        pass
                    else:
                        md = str("S" + str(r) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                        g.write(md)

            with open(nmoc_metadata, 'a') as g:
                if ((len(x[8]) > 1) or ('-' in x[3] and '_(' not in x[3])):
                    pass
                else:
                    if z == 1:
                        g.write(meta_first_line)
                    if (d.match(x[3])):
                        pass
                    else:
                        md = str("S" + str(z) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                        g.write(md)

            with open(wt_vert_metadata, 'a') as g:
                if (x[6] != "Chordata") or d.match(x[3]) or ((len(x[8]) > 1) or ('-' in x[3] and '_(' not in x[3])):
                    pass
                else:
                    if wt_vert == 1:
                        g.write(meta_first_line)

                    md = str("S" + str(wt_vert) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                    g.write(md)

            with open(mut_metadata, 'a') as g:
                if ((len(x[8]) > 1) or ('-' in x[3] and '_(' not in x[3])):
                    if mut == 1:
                        g.write(meta_first_line)
                    if (d.match(x[3])):
                        pass
                    else:
                        md = str("M" + str(mut) + "\t" + str(lmax).strip()) + "\t" + str(x[0]).strip().replace(' ', '') + "_" + str(x[1]).strip().replace(' ', '') + "\t" + str(x[2]).strip() + "\t" + str(x[6]).strip() + "\t" + str(x[7]).strip().replace(' ', '') + '\t' + x[3].strip() + "\t" + str(x[8]).strip() + "\t" + str(x[5]).strip().replace(' ', '') + "\t" + str(x[9]).strip() + "\n"
                        g.write(md)
                else:
                    pass


    # Calculate average lmax and update for sequences flagged as haivng redundant entries with conflicting lmax values in the database
    for dsplit in meta_data_list:
        file = fileinput.input(dsplit, inplace=True)
        for line in file:
            modify = False
            fields = line.strip().split("\t")
            lmax = fields[1]
            accession = fields[6]
            for key in conflicting_lmax.keys():
                if accession == key:
                    modify = True
            if modify:
                conflicting_lmax[accession]
                avg_lmax = float(sum(conflicting_lmax[accession])) / float(len(conflicting_lmax[accession]))
                fields[1] = str(avg_lmax)
            print("\t".join(fields))  # Print the modified line, which overwrites the original line in the file
        file.close()
    # Creates an additional sequence file where all the mutants are
    # added to the bottom of the wild-type sequences so they can be aligned for WT model test.
    mut_only = open(mut_output).readlines()
    mut_nmoc = f'{seq_report_dir}/wt_mut_added.txt'
    shutil.copy(nmoc_output, mut_nmoc)
    x = 0
    for lines in mut_only:
        if x <= 1:
            if x == 0:
                with open(mut_nmoc, 'a') as m:
                    m.write('\n')
            else:
                pass
            x += 1
        else:
            with open(mut_nmoc, 'a') as m:
                m.write(lines)

    return seq_report_dir, data_split_list, meta_data_list, mnm_meta_list, mnm_meta_shorthand

# seq_report_dir, data_split_list, meta_data_list, mnm_meta_list, mnm_meta_shorthand = extract_vpod_datasets(mydb)
# print(f"VPOD datasets extracted to: {seq_report_dir}")

from deepBreaks.preprocessing import read_data
import pandas as pd
import shutil

def get_mnm_datasets(seq_report_dir, mnm_data, mnm_meta_list, mnm_meta_shorthand, meta_data_list, data_split_list):
    for meta, short_hand in zip(mnm_meta_list,mnm_meta_shorthand):
        
        meta_df = read_data(meta, seq_type = None, is_main=False)
        mnm_data_copy = mnm_data.copy()
        mnm_data_copy.drop_duplicates(subset=['Protein'], keep='first', inplace=True)
        if short_hand == 'vert'  or short_hand == 'wt_vert':
            mnm_data_copy = mnm_data_copy[mnm_data_copy['Phylum'] == 'Chordata']
        elif short_hand == 'inv':
            mnm_data_copy = mnm_data_copy[mnm_data_copy['Phylum'] != 'Chordata']
        else:
            pass
            
        mnm_data_copy = mnm_data_copy[~mnm_data_copy['Protein'].isin(meta_df['Protein'].to_list())]
        mnm_ids = list(range(len(meta_df),len(meta_df)+len(mnm_data_copy['Protein'].to_list())))
        
        combined_ids = meta_df.index.to_list()
        for ids in mnm_ids:
            combined_ids.append(f'MNM{ids}')

        combined_prot = meta_df['Protein'].to_list()
        for prot in mnm_data_copy['Protein']:
            combined_prot.append(prot)
            
        combined_lmax = meta_df['Lambda_Max'].to_list()
        for lmax in mnm_data_copy['LambdaMax']:
            combined_lmax.append(lmax)

        combined_acc = meta_df['Accession'].to_list()
        for acc in mnm_data_copy['Accession']:
            combined_acc.append(acc)

        combined_sp = meta_df['Species'].to_list()
        mnm_sp_list = mnm_data_copy['Genus'] + '_' + mnm_data_copy['Species']
        for sp in mnm_sp_list:
            combined_sp.append(sp)
            
        merged_df = pd.DataFrame(columns=['Accession','Species','Lambda_Max', 'Protein'], index=combined_ids)
        merged_df['Accession'] = combined_acc
        merged_df['Species'] = combined_sp
        merged_df['Lambda_Max']= combined_lmax
        merged_df['Protein']= combined_prot
        merged_df_file = f'{seq_report_dir}/{short_hand}_mnm_meta.csv'
        merged_df.to_csv(path_or_buf=f'./{merged_df_file}', index=True)
        
        vpod_seqs = f'{seq_report_dir}/{short_hand}.txt'
        vpod_mnm_seqs = f'{seq_report_dir}/{short_hand}_mnm.txt'
        shutil.copy(vpod_seqs , vpod_mnm_seqs)
        mnm_seqs = mnm_data_copy['Protein'].to_list()

        with open(vpod_mnm_seqs, 'a') as f:
            for id, seq in zip(mnm_ids, mnm_seqs):
                f.write(f'>MNM{id}\n{seq}\n')
                
        data_split_list.append(vpod_mnm_seqs)
        meta_data_list.append(merged_df_file)
        
    return data_split_list, meta_data_list


from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import os  # Import os module if not already imported

def perform_mafft_alignment(data_split_list, mafft_exe_path):
    """
    Performs multiple sequence alignment using MAFFT on a list of input files.

    Args:
        data_split_list (list): A list of paths to the input sequence files (e.g., FASTA files).
        mafft_exe_path (str): The path to the MAFFT executable.

    Returns:
        list: A list of paths to the aligned output files ready for input into deepBreaks.
    """
    mafft_output_list = []
    for data in data_split_list:
        output = f'./{data.split(".")[0]}.{data.split(".")[1]}_aligned.txt'
        mafft_cline = MafftCommandline(mafft_exe_path, input=f'./{data}')
        #print(mafft_cline)
        stdout, stderr = mafft_cline()

        with open(output, "w") as handle:
            handle.write(stdout)
            #print(handle)
        align = AlignIO.read(output, "fasta")
        mafft_output_list.append(f'{output}')
        
    deep_breaks_input_data = []
    for input_file in mafft_output_list:
        #print(item)
        output_file = f'.{input_file.split(".")[1]}.{input_file.split(".")[2]}_VPOD_1.2_het.fasta'
        deep_breaks_input_data.append(output_file)
        with open(output_file, 'w') as file:
                with open(input_file, 'r') as infile: # Use with open for file handling
                    lines = infile.readlines()
                    m = 0
                    for k, line in enumerate(lines): #Enumerate lines for index tracking
                        snip = line.strip() # Strip whitespace for cleaner processing
                        if '>' in snip:
                            if m > 0: # Check if it is not the first sequence header
                                file.write("\n") # Add newline before new header (except for first one)
                            m += 1
                            file.write(snip + "\n") # Write header and newline immediately
                        else:
                            file.write(snip) # Write sequence line without adding newline - original code behaviour

    print(deep_breaks_input_data)
    return deep_breaks_input_data
