#!/usr/bin/python3
import argparse
import logging
import os
import shutil
import pathlib
import sys
import string
import pandas
import random
import io
import re
import math
from io import BytesIO
from subprocess import Popen, PIPE
import datetime
import multiprocessing as mp
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import decouple
import numpy as np




script_path = os.getcwd()
config = decouple.Config(decouple.RepositoryEnv(script_path + '/conf.env'))
SignalP4_path = config('SIGNALP4')
SignalP5_path = config('SIGNALP5')
Phobius_path = config('PHOBIUS')
LipoP_path = config('LIPOP')
TMHMM_path = config('TMHMM')
TatP_path = config('TATP')
Interpro_path = config('INTERPRO')
PERL5LIB = config('PERL5LIB')

#####

__author__ = 'StefanoGrasso'

parser = argparse.ArgumentParser(description='Script to perform a proteomic based consensus prediction of '
                                             'protein localization in Gram+. Written by Stefano Grasso (c).',
                                 prog='GP4')
parser.add_argument('-i', '--input', help='Input file', required=True)
parser.add_argument('-o', '--output', help='Output folder name', required=False, default='Results')
parser.add_argument('-n', '--name', help='Output file name', required=False)
parser.add_argument('-t', '--threads', help='Number of threads to be used. Default: all.', required=False)
parser.add_argument('-r', '--run', help='\'data\' to generate data; \'pred\' requires generated data as input; \'all\' (default)',
                    choices=['data', 'pred', 'all'] ,required=False, default='all')
parser.add_argument('-v', '--version', action='version', version='%(prog)s v. 1.0-gamma')
args = parser.parse_args()

t_start = datetime.datetime.now()




if args.run == 'data' or args.run == 'all':
    full_path = os.path.abspath((args.input))
    filename_w_ext = os.path.basename(full_path)
    filename = os.path.splitext(filename_w_ext)[0]
    output_path = str(os.path.abspath((args.output))) + '/'
    if args.name:
        log_name = str(args.name)
        filename = str(args.name)
    else:
        log_name=filename

    logging.basicConfig(filename='GP4_'+log_name+'.log', level=logging.DEBUG, filemode='w')
    logging.info("Script sarted at " + str(t_start))

    if os.path.isfile(SignalP4_path):
        logging.info("SignalP4.0 executable correct")
    else:
        logging.error("SignalP4.0 executable not found"), sys.exit('SignalP4.0 executable not found')
    if os.path.isfile(SignalP5_path):
        logging.info("SignalP5.0 executable correct")
    else:
        logging.error("SignalP5.0 executable not found"), sys.exit('SignalP5.0 executable not found')
    if os.path.isfile(Phobius_path):
         logging.info("Phobius executable correct")
    else: logging.error("Phobius executable not found"), sys.exit('Phobius executable not found')
    if os.path.isfile(LipoP_path):
        logging.info("LipoP executable correct")
    else:
        logging.error("LipoP executable not found"), sys.exit('LipoP executable not found')
    if os.path.isfile(TMHMM_path):
        logging.info("TMHMM executable correct")
    else:
        logging.error("TMHMM executable not found"), sys.exit('TMHMM executable not found')
    if os.path.isfile(TatP_path):
        logging.info("TatP executable correct")
    else:
        logging.error("TatP executable not found"), sys.exit('TatP executable not found')
    if os.path.isfile(Interpro_path):
        logging.info("Interpro executable correct")
    else:
        logging.error("Interpro executable not found"), sys.exit('Interpro executable not found')
    if os.path.isfile(args.input):
        logging.info("Input file exists")
    else:
        logging.error("Input file not found"), sys.exit('Input file not found')
    logging.info("Input file is: " + full_path + ". File name is: " + filename)
    logging.info("Output will be in: " + output_path)
    logging.info("Number of CPUs available: " + str(os.cpu_count()))
    if args.threads:
        try:
            n_thread = int(args.threads)
        except:
            logging.error("Invalid number of CPUs: value inserted not a number. ")
            exit()
        if int(args.threads) < os.cpu_count():
            ncpu = int(args.threads)
            logging.info("Number of CPUs to be used: " + str(ncpu))
        else:
            logging.error("Invalid number of CPUs: too many. ")
            exit()
    else:
        ncpu = int(os.cpu_count())

    try:
        pathlib.Path(output_path + filename).mkdir(parents=True, exist_ok=False)
        new_foldername = filename
    except:
        logging.warning("Folder already exists, new name is being given")
        rnd_num = ''.join(random.choices(string.ascii_letters + string.digits, k=5))
        new_foldername = filename + '_' + rnd_num
        logging.info("New folder and output file name is: " + new_foldername)
        pathlib.Path(output_path + new_foldername).mkdir(parents=True, exist_ok=False)

    logging.info("Generating new fasta with shorter id")
    shortid_file = output_path + new_foldername + '/' + new_foldername + '.shortid.txt'


    old_ids = []
    with open(args.input, 'r') as handle:
        with open(shortid_file, 'w') as h:
            newseqs = []
            for record in SeqIO.parse(handle, 'fasta', alphabet=(IUPAC.extended_protein)):
                old_ids.append(record.id)
                newrecord = SeqIO.SeqRecord(record.seq, record.id[-19:], '', '')
                newseqs.append(newrecord)
            SeqIO.write(newseqs, h, 'fasta')
    logging.info("New fasta saved to " + shortid_file)
    with open(output_path + new_foldername + '/ID_' + new_foldername + '.id', 'w') as handle:
        for id in old_ids:
            handle.write("%s\n" % id)


elif args.run == 'pred':
    full_path = os.path.abspath((args.input))
    output_path = full_path.rsplit('/',1)[0] +'/'
    new_foldername = os.path.basename(full_path)
    logging.basicConfig(filename='GP4_'+new_foldername+'.log', level=logging.DEBUG, filemode='w')
    log_name=new_foldername
    logging.info("Script sarted at " + str(t_start))
    if os.path.isdir(full_path):
        logging.info("Input folder is: " + new_foldername + ". Full path is: " + full_path)
        shortid_file = full_path + new_foldername + '/' + new_foldername + '.shortid.txt'

    else:
        logging.error("Invalid input folder. Input is not a folder or does not exist.")
        exit()
else:
    logging.error("--run option not allowed.")
    exit()

print(output_path + new_foldername + '/SignalP4_' + new_foldername + '.csv')
outfile_SigP =  output_path + new_foldername + '/SignalP4_' + new_foldername + '.csv'
outfile_SigP5 = output_path+ new_foldername + '/SignalP5_' + new_foldername + '.log'
results_sigP5 = output_path + new_foldername + '/' + new_foldername + '_summary.signalp5'
outfile_Pho = output_path + new_foldername + '/Phobius_' + new_foldername + '.csv'
outfile_LipoP = output_path + new_foldername + '/LipoP_' + new_foldername + '.csv'
outfile_TMHMM = output_path + new_foldername + '/TMHMM_' + new_foldername + '.csv'
outfile_PCB = output_path + new_foldername + '/ProtCompB_' + new_foldername + '.csv'
outfile_Predisi = output_path + new_foldername + '/PrediSi_' + new_foldername + '.csv'
outfile_TatP = output_path + new_foldername + '/TatP_' + new_foldername + '.csv'
outfolder_Interpro = output_path + new_foldername + '/Interpro_' + new_foldername
outlog_Interpro = outfolder_Interpro + '/Log_' + new_foldername + '.log'
outfile_results = output_path + new_foldername + '/Final_' + new_foldername + '.csv'


def strip(text):
    try:
        return text.strip()
    except:
        return text

def get_val(text):
    try:
        return text.split('=')[1]
    except:
        return text


def SIGNALP(path):
    logging.info("Started SignalP4.1")
    SignalP4_process = Popen([SignalP4_path, '-t', 'gram+', '-f', 'short', path], stdout=PIPE, stderr=PIPE)
    stdout, stderr = SignalP4_process.communicate()
    if stderr: logging.warning(stderr)
    reader = io.BufferedReader(io.BytesIO(stdout))
    wrapper = io.TextIOWrapper(reader)
    lines = wrapper.read().splitlines()[1:]
    stdout = [x.split() for x in lines]
    headers = stdout.pop(0)[1:]
    df = pandas.DataFrame(stdout, columns=headers)
    with open(outfile_SigP, 'w') as handle:
        df.to_csv(handle, sep='\t')
    logging.info("SignalP4.1 results saved to " + outfile_SigP)


def PHOBIUS(path):
    logging.info("Started Phobius")
    Phobius_process = Popen([Phobius_path, '-short', path], stdout=PIPE, stderr=PIPE)
    stdout, stderr = Phobius_process.communicate()
    if stderr: logging.warning(stderr)
    reader = io.BufferedReader(io.BytesIO(stdout))
    wrapper = io.TextIOWrapper(reader)
    lines = wrapper.read().splitlines()
    stdout = [x.split() for x in lines]
    headers = stdout.pop(0)
    headers[0:2] = ['_'.join(headers[0:2])]
    df = pandas.DataFrame(stdout, columns=headers)
    with open(outfile_Pho, 'w') as handle:
        df.to_csv(handle, sep='\t')
    logging.info("Phobius results saved to " + outfile_Pho)


def LIPOP(path):
    logging.info("Started LipoP")
    LipoP_process = Popen([LipoP_path, '-short', path], stdout=PIPE, stderr=PIPE)
    stdout, stderr = LipoP_process.communicate()
    if stderr: logging.warning(stderr)
    reader = io.BufferedReader(io.BytesIO(stdout))
    wrapper = io.TextIOWrapper(reader)
    lines = wrapper.read().splitlines()
    stdout = [x.split() for x in lines]
    headers = ['TODROP', 'ID', 'Localization', 'Score', 'Margin', 'Cl_Site', 'Pos+2']
    df = pandas.DataFrame(stdout)
    df.columns = headers[:(len(df.columns))]
    df = df.drop(['TODROP'], axis=1)
    with open(outfile_LipoP, 'w') as handle:
        df.to_csv(handle, sep='\t')
    logging.info("LipoP results saved to " + outfile_LipoP)


def TMHMM(path):
    logging.info("Started TMHMM")
    TMHMM_process = Popen([TMHMM_path, '-short', '-noplot', path], stdout=PIPE, stderr=PIPE)
    stdout, stderr = TMHMM_process.communicate()
    if stderr: logging.warning(stderr)
    stdout_f = BytesIO(stdout)
    df = pandas.read_csv(stdout_f, sep=';', skiprows=0, names=['ID', 'Len', 'AA_TM', 'First60', 'Num_Hel', 'Topology'])
    with open(outfile_TMHMM, 'w') as handle:
        df.to_csv(handle, sep='\t')
    logging.info("TMHMM results saved to " + outfile_TMHMM)


def PREDISI(path):
    logging.info("Started Predisi")
    Predisi_process = Popen(['curl', '--form', 'Input=select', '--form', 'Datei=@'+path,
                             '--form', 'matrix=Gram-positive', '--form', 'output=CSV', '--form', 'delimiter=,',
                             'http://www.predisi.de/predisi/startprediction'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = Predisi_process.communicate()
    if stderr: logging.warning(stderr)
    stdout_f = BytesIO(stdout)
    df = pandas.read_csv(stdout_f, sep=',', skiprows=4, index_col=None)
    with open(outfile_Predisi, 'w') as handle:
        df.to_csv(handle, sep='\t')
    logging.info("Predisi results saved to" + outfile_Predisi)


def PROTCOMPB(path):
    logging.info("Started ProtCompB")
    PCB_process = Popen(['curl', '--form', 'FILE=@' + path, '--form', 'gramm=1',
                         'http://linux1.softberry.com/cgi-bin/programs/proloc/pcompb.pl'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = PCB_process.communicate()
    if stderr: logging.warning(stderr)
    reader = io.BufferedReader(io.BytesIO(stdout))
    wrapper = io.TextIOWrapper(reader)
    lines = wrapper.read().splitlines()
    seq_ids = [s for s in lines if 'Seq name' in s]
    seqs = []
    for l in seq_ids:
        seqs.append([re.split(r' |,', l)[2]])
    results = [s for s in lines if 'Integral Prediction' in s]
    for i, l in zip(seqs, results):
        i.append((re.split(' ', l)[5]))
        i.append((re.split(' ', l)[-1]))
    df = pandas.DataFrame(seqs, columns=['ID', 'Prediction', 'Score'])
    with open(outfile_PCB, 'w') as handle:
        df.to_csv(handle, sep='\t')
    logging.info("ProtCompB results saved to " + outfile_PCB)


def TATP(path):
    logging.info("Started TatP")
    TatP_process = Popen([TatP_path, '-short', '-trunc', '100', path], stdout=PIPE, stderr=PIPE,
                         env={'PERL5LIB': PERL5LIB})
    stdout, stderr = TatP_process.communicate()
    if stderr: logging.warning(stderr)
    reader = io.BufferedReader(io.BytesIO(stdout))
    wrapper = io.TextIOWrapper(reader)
    lines = wrapper.read().splitlines()
    seq_ids = [[s, lines.index(s)] for s in lines if '>' in s]
    seqs = []
    j = 0
    for l, i in seq_ids:
        seqs.append([str(re.split(r' ', l)[0])[1:]])
        if "Most likely cleavage site" in lines[i + 7]:
            cl_b1, cl_a1, consensus = re.split(" ", lines[i + 7])[7], re.split(" ", lines[i + 7])[9][:-1], \
                                      re.split(" ", lines[i + 7])[-1]
            if "as Tat motif starting at" in lines[i + 8]:
                tat_motif, tat_pos = re.split(" ", lines[i + 8])[2], re.split(" ", lines[i + 8])[-2]
            elif "Potential Tat signal peptide but No Tat motif" in lines[i + 8]:
                tat_motif, tat_pos = 'No Tat motif', 'No Tat motif'
            else:
                tat_motif, tat_pos = 'No Tat motif', 'No Tat motif'
            seqs[j].extend(['Y', cl_b1, cl_a1, consensus, tat_motif, tat_pos])
        elif "Used regex" in lines[i + 7]:
            seqs[j].extend(['N', '', '', '', '', ''])
        else:
            seqs[j].extend(['ERROR', '', '', '', '', ''])
        j += 1
    df = pandas.DataFrame(seqs, columns=['ID', 'SignalPeptide', 'Cl.Site-1', 'Cl.Site+1', 'Cl.Seq', 'TatPMotif',
                                         'TatPos'])
    with open(outfile_TatP, 'w') as handle:
        df.to_csv(handle, sep='\t')
    logging.info("TatP results saved to " + outfile_TatP)


def INTERPRO(path):
    logging.info("Started InterproScan")
    pathlib.Path(outfolder_Interpro).mkdir(parents=True, exist_ok=False)
    Interpro_process = Popen([Interpro_path, '-d ', outfolder_Interpro, '-iprlookup', '-pa', '-etra', '-goterms', '-i', path], stdout=PIPE,
        stderr=PIPE, env={'PERL5LIB': PERL5LIB})
    stdout, stderr = Interpro_process.communicate()
    with open(outlog_Interpro, 'wb') as handle:
        handle.write(stdout)
    if stderr: logging.warning(stderr)
    logging.info("InterproScan results saved to " + outfolder_Interpro)


def Interpro_arrange(row, list_domains, df, type):
    if type == 'IPR':
        annot_col = 'IPR_annotation_Interpro'
        descr_col = 'IPR_description_Interpro'
    elif type == 'NI':
        annot_col = 'Signature Accession_Interpro'
        descr_col = 'Signature Description_Interpro'
    elif type == 'GO':
        annot_col = 'GO_annotation_Interpro'
    for index_ipr, row_ipr in df.iterrows():
        if row['name_SigP4'] == row_ipr['ID_Interpro']:
            list_idx = pandas.DataFrame(columns=['ID','DOMAIN','IPR','DESCRIPTION'])
            for i, s in enumerate(row_ipr[annot_col]):
                if not is_nan(s):
                    dom, ipr = searchDict(list_domains, str(s))
                    if not dom == False:
                        if type == 'IPR' or type == 'NI':
                            line = {'ID': row_ipr['ID_Interpro'],'DOMAIN': dom, 'IPR': s, 'DESCRIPTION': row_ipr[descr_col][i]}
                        elif type == 'GO':
                            line = {'ID': row_ipr['ID_Interpro'],'DOMAIN': dom, 'IPR': s}
                        list_idx = list_idx.append(line, ignore_index=True)
            if list_idx.empty:
                if type == 'IPR' or type == 'NI':
                    return np.nan, np.nan, np.nan
                elif type == 'GO':
                    return np.nan, np.nan
            else:
                if type == 'IPR' or type == 'NI':
                    return list_idx['DOMAIN'].tolist(), list_idx['IPR'].tolist(), list_idx['DESCRIPTION'].tolist()
                elif type == 'GO':
                    return list_idx['DOMAIN'].tolist(), list_idx['IPR'].tolist()

def Predict(row):
    sec = 0
    tm = 0
    lipo = 0
    cyt = 0
    sp_cl = []
    lipo_cl = []
    surf_b = False
    lant_b =False
    lipo_b =False
    cw_mp = False
    cw_nc = False
    cw_spore = False
    tat = False
    sec_b =False
    pilin_b = False
    tm_b = False
    cyt_b = False

#SignalP4.1
    if row['?_SigP4'] == 'Y':
        sec +=1
        sp_cl.append(int(row['pos.2_SigP4']))
    elif row['?_SigP4'] == 'N':
        cyt +=1
#SignalP5
    if row['Prediction_SigP5'] == 'SP(Sec/SPI)' or row['Prediction_SigP5'] == 'TAT(Tat/SPI)':
        sec +=1
        cl = re.search(r':(.*?)\.', str(row['CS Position_SigP5'])).group(1)
        cl = cl.split('-', 1)[1]
        sp_cl.append(int(cl.strip()))
        sp_cl.append(int(cl.strip())) ### duplicated on purpose. if SigP5 then it has double weight.
    elif row['Prediction_SigP5'] == 'LIPO(Sec/SPII)':
        lipo +=1
        cl = re.search(r':(.*?)\.', str(row['CS Position_SigP5'])).group(1)
        cl = cl.split('-', 1)[1]
        lipo_cl.append(int(cl.strip()))
    elif row['Prediction_SigP5'] == 'OTHER':
        cyt +=1
        tm +=1
#Phobius
    if row['SP_Phobius'] == 'Y':
        sec +=1
        cl = re.search(r'/(.*?)o', str(row['PREDICTION_Phobius'])).group(1)
        sp_cl.append(int(cl.strip()))
    elif row['SP_Phobius'] == '0':
        cyt +=1
    if int(row['TM_Phobius']) < 3:
        sec += 1
        tm += 1
    elif int(row['TM_Phobius']) == 0:
        cyt += 1
    elif int(row['TM_Phobius']) >= 3:
        tm += 2
#Predisi
    if row['Signal Peptide ?_PrediSi'] == 'Y':
        sec +=1
        cl = int(row['Cleavage Position_PrediSi']) + 1
        sp_cl.append(cl)
    elif row['Signal Peptide ?_PrediSi'] == 'N':
        cyt +=1
#TatP
    if row['SignalPeptide_TatP'] == 'Y':
        sec += 1
        cl = int(float(str(row['Cl.Site+1_TatP']).strip()))
        sp_cl.append(cl)
    if not row['TatPMotif_TatP'] == 'No Tat motif' and not is_nan(row['TatPMotif_TatP']):
        sec += 1
        tat = True
#LipoP
    if row['Localization_LipoP'] == 'SpI':
        sec += 1
        cl = int(float(str(row['Cl_Site_LipoP']).split('-')[1]))
        sp_cl.append(cl)
    elif row['Localization_LipoP'] == 'SpII':
        lipo += 1
        cl = int(str(row['Cl_Site_LipoP']).split('-')[1])
        lipo_cl.append(cl)
    elif row['Localization_LipoP'] == 'CYT':
        cyt += 1
    elif row['Localization_LipoP'] == 'TMH':
        tm += 1
#TMHMM
    if int(row['Num_Hel_TMHMM']) == 0:
        cyt += 1
    elif int(row['Num_Hel_TMHMM']) < 3:
        sec += 1
        tm += 1
    elif int(row['Num_Hel_TMHMM']) >= 3:
        tm += 2
#PCB
    if row['Prediction_PCB'] == 'Secreted':
        sec += 1
    elif row['Prediction_PCB'] == 'Cytoplasmic':
        cyt += 1
    elif row['Prediction_PCB'] == 'Membrane':
        tm += 1
#IPR
    ni_dom = row['Domain_type_ni']
    if not is_nan(ni_dom):
        for dom in ni_dom:
            if dom == 'CW_mp_ni':
                cw_mp = True
                cyt -= 1
    ipr_dom = row['Domain_type']
    if not is_nan(ipr_dom):
        for dom in ipr_dom:
            if dom == 'CW_mp':
                cw_mp = True
                cyt -= 1
            elif dom == 'CW_nc':
                cw_nc = True
                cyt -= 1
            elif dom == 'CW_spore':
                cw_spore = True
                cyt -= 1
            elif dom == 'Lant':
                sec += 1
                cyt -= 1
                lant_b = True
            elif dom == 'Lipo':
                lipo += 1
                lipo_b = True
            elif dom == 'Pilin':
                sec += 1
                cyt -= 1
                pilin_b= True
            elif dom == 'Sec':
                sec += 1
                cyt -= 1
                sec_b = True
            elif dom == 'Surf':
                sec += 1
                cyt -= 1
                surf_b = True
            elif dom == 'Tat':
                sec += 1
                tat = True

    go_dom = row['GO_annotation_Interpro']
    if not is_nan(go_dom):
        for dom in go_dom:
            if dom == 'CW':
                cw_nc = True
                cyt -= 1
            elif dom == 'CW_spore':
                sec += 1
                cyt -= 1
                cw_spore = True
            elif dom == 'TM':
                tm += 1
                tm_b = True
            elif dom == 'CYTO':
                cyt += 1
                cyt_b = True
            elif dom == 'Pilin':
                sec += 1
                cyt -= 1
                pilin_b= True
            elif dom == 'Sec':
                sec += 1
                cyt -= 1
                sec_b = True
            elif dom == 'Surf':
                sec += 1
                cyt -= 1
                surf_b = True
#Returns
    result = []
    note = []
    if sec_b or lant_b or sec >= 3:
        result.append('EXTRA')
        if len(sp_cl) >= 2:
            cl_site = most_common(sp_cl)
            note.append('Contains a SP (Sec/SpI), most likely first a.a. of the mature protein: ' + str(cl_site) + '.')
        if tat == True:
            note.append('Secreted via TAT, motif ' + str(row['TatPMotif_TatP']) + '.')
        if sec_b == True:
            note.append('Secreted via a non-canonical secretion pathway.')
        if lant_b == True:
            note.append('Short secreted peptide, could be a bacteriocin, lantibiotic, or similar.')
    if cw_mp or cw_nc or cw_spore or surf_b:
        result.append('CW')
        cyt += -3
        if cw_mp:
            note.append('Covalently attached to the CW.')
        if cw_nc:
            note.append('Probably interacts, not covalently, with the CW.')
        if cw_spore:
            note.append('Should localize at division/spore septum or spore cortex or coat.')
        if surf_b == True:
            note.append('Secreted and displayed on the surface, possibly an S-layer protein.')
    if tm_b or tm >= 4:
        result.append('TM')
    if lipo_b or lipo >= 2:
        result.append('LIPO')
        if len(lipo_cl) >= 1:
            cl_site = most_common(lipo_cl)
            for n in note:
                if '(Sec/SpI)' in n:
                    note.remove(n)
            note.append('Contains a lipoprotein SP (Sec/SpII), most likely first a.a. of the mature protein: ' + str(cl_site) + '.')
    if cyt_b or cyt >= 5:
        if cyt_b:
            result = ['CYTO']
        elif len(result) == 0:
            result.append('CYTO')
    if pilin_b == True:
        note.append('Part of pilin, flagellum, fimbriae, competence protein or similar.')

    if len(result) > 0:
        if len(result) >= 2:
            if (('EXTRA' in result) and ('CW' in result)) or (('EXTRA' in result) and ('LIPO' in result)):
                result.remove('EXTRA')
            if (('CYTO' in result) and ('LIPO' in result)) or (('CYTO' in result) and ('TM' in result)):
                result.remove('CYTO')
        if len(result) == 2:
            if ('EXTRA' in result) and ('CYTO' in result):
                result = ['Unknown']
            if ('EXTRA' in result) and ('TM' in result):
                result = ['Unknown']
        elif len(result) == 3:
            if ('EXTRA' in result) and ('TM' in result) and ('CYTO' in result):
                result = ['Unknown']
        return ' '.join(result), ' '.join(note)
    else: return 'Unknown', ''


def SIGNALP5(path):
    logging.info("Started SignalP5")
    SigP5_process = Popen([SignalP5_path, '-fasta', path, '-org', 'gram+', '-prefix', output_path + new_foldername+'/'+new_foldername], stdout=PIPE, stderr=PIPE)
    stdout, stderr = SigP5_process.communicate()
    if stderr: logging.warning(stderr)
    with open(outfile_SigP5, 'wb') as handle:
        handle.write(stdout)
    logging.info("SignalP5 log saved to " + outfile_SigP5)


def LOADLIST(infile):
    domains = []
    with open(infile) as handle:
        for line in handle:
            domains.append(line.strip())
    domains = list(set(domains))
    return domains

def searchDict(myDict, lookup):
    for key, value in myDict.items():
        for v in value:
            if str(lookup) in str(v):
                return key, v
    return False, False

def is_nan(val):
    try:
        return math.isnan(float(val))
    except:
        return False

def most_common(lst):
    return max(set(lst), key=lst.count)

def split_GO(row):
    go_list = []
    for i, s in enumerate(row['GO_annotation_Interpro']):
        if not is_nan(s):
            if 'GO' in s:
                try:
                    go = str(s).strip().split('|')
                    go_list.extend(go)
                except:
                    go_list.extend(s)
    if go_list:
        return go_list
    else: return [np.nan]


if args.run == 'data' or args.run == 'all':

    # if args.threads:
    #     if type(args.threads) == int and args.threads < os.cpu_count():
    #         ncpu = args.threads
    #     else:
    #         ncpu = os.cpu_count()
    # else:
    #     ncpu = os.cpu_count()

    pool = mp.Pool(processes=ncpu)

    SigP5_out = pool.apply_async(SIGNALP5, [shortid_file])
    SigP_out = pool.apply_async(SIGNALP, [shortid_file])
    Pho_out = pool.apply_async(PHOBIUS, [shortid_file])
    Predisi_out = pool.apply_async(PREDISI, [shortid_file])
    LipoP_out = pool.apply_async(LIPOP, [shortid_file])
    TMHMM_out = pool.apply_async(TMHMM, [shortid_file])
    PCB_out = pool.apply_async(PROTCOMPB, [shortid_file])
    TatP_out = pool.apply_async(TATP, [shortid_file])
    Interpro_out = pool.apply_async(INTERPRO, [shortid_file])

    pool.close()
    pool.join()

    logging.info("All data are generated.")

if args.run == 'pred' or args.run == 'all':

    logging.info("Now final predictions are being done.")

    with open(outfile_SigP, 'r') as handle:
        results_df = pandas.read_csv(handle, sep='\t', index_col=0, converters={'name': strip})
        results_df.columns = results_df.columns.map(lambda x: str(x) + '_SigP4')

    with open(results_sigP5, 'r') as handle:
        SigP5_df = pandas.read_csv(handle, sep='\t', skiprows=1, index_col=None, converters={'# ID': strip})
        SigP5_df.columns = SigP5_df.columns.map(lambda x: str(x) + '_SigP5')
        if len(results_df.index) == len(SigP5_df.index):
            results_df = pandas.concat([results_df, SigP5_df], axis=1 )  #left_on='name_SigP4', right_on='# ID_SigP5'
        else:
            logging.error('SigP5 ERROR, protein num does not match')

    with open(outfile_Pho) as handle:
        Pho_df = pandas.read_csv(handle, sep='\t', index_col=0, converters={'SEQENCE_ID': strip})
        Pho_df.columns = Pho_df.columns.map(lambda x: str(x) + '_Phobius')
        results_df = results_df.merge(Pho_df, how='left', left_on='name_SigP4', right_on='SEQENCE_ID_Phobius')

    with open(outfile_Predisi) as handle:
        Predisi_df = pandas.read_csv(handle, sep='\t', index_col=0, converters={'FASTA-ID': strip})
        Predisi_df.columns = Predisi_df.columns.map(lambda x: str(x) + '_PrediSi')
        results_df = results_df.merge(Predisi_df, how='left', left_on='name_SigP4', right_on='FASTA-ID_PrediSi')

    with open(outfile_TatP) as handle:
        Tat_df = pandas.read_csv(handle, sep='\t', index_col=0, converters={'ID': strip})
        Tat_df.columns = Tat_df.columns.map(lambda x: str(x) + '_TatP')
        results_df = results_df.merge(Tat_df, how='left', left_on='name_SigP4', right_on='ID_TatP')

    with open(outfile_LipoP) as handle:
        Lipo_df = pandas.read_csv(handle, sep='\t', index_col=0, converters={'ID': strip, 'Score': get_val, 'Margin': get_val, 'Cl_Site': get_val})
        Lipo_df.columns = Lipo_df.columns.map(lambda x: str(x) + '_LipoP')
        results_df = results_df.merge(Lipo_df, how='left', left_on='name_SigP4', right_on='ID_LipoP')

    with open(outfile_TMHMM) as handle:
        TMHMM_df = pandas.read_csv(handle, sep='\t', index_col=0, converters={'ID': strip, 'Len': get_val, 'AA_TM': get_val,
                                                                              'First60': get_val, 'Num_Hel': get_val, 'Topology': get_val})
        TMHMM_df.columns = TMHMM_df.columns.map(lambda x: str(x) + '_TMHMM')
        results_df = results_df.merge(TMHMM_df, how='left', left_on='name_SigP4', right_on='ID_TMHMM')

    with open(outfile_PCB) as handle:
        PCB_df = pandas.read_csv(handle, sep='\t', index_col=0, converters={'SeqID': strip})
        PCB_df.columns = PCB_df.columns.map(lambda x: str(x) + '_PCB')
        PCB_df['ID_PCB'] = PCB_df['ID_PCB'].astype(str)
        results_df = results_df.merge(PCB_df, how='left', left_on='name_SigP4', right_on='ID_PCB')

    Domains = {'Pilin': LOADLIST('Domains/Pilin.ipr'),
    'CW_mp': LOADLIST('Domains/CW_mp.ipr'),
    'CW_nc': LOADLIST('Domains/CW_nc.ipr'),
    'CW_spore': LOADLIST('Domains/CW_spore.ipr'),
    'Lant': LOADLIST('Domains/Lant.ipr'),
    'Lipo': LOADLIST('Domains/Lipo.ipr'),
    'Sec': LOADLIST('Domains/Sec.ipr'),
    'Surf': LOADLIST('Domains/Surf.ipr'),
    'Tat': LOADLIST('Domains/Tat.ipr')}

    Domains_ni = {'CW_mp_ni': LOADLIST('Domains/CW_mp.iprni')}

    GO_dict = {'Sec': ['GO:0005576'], 'CW': ['GO:0005618'], 'CW_spore': ['GO:0031160'], 'TM': ['GO:0005886', 'GO:0016021'],
               'CYTO': ['GO:0005829', 'GO:0005737', 'GO:0009295','GO:0043590'], 'Surf': ['GO:0009986', 'GO:0009897'], 'Pilin': ['GO:0009289', 'GO:0009288']}


    outfile_Interpro = outfolder_Interpro + '/' + os.path.basename(shortid_file) + '.tsv'
    headers = ['ID', 'MD5 digest', 'Length', 'Analysis', 'Signature Accession', 'Signature Description',
               'Start location',
               'Stop location', 'Score', 'Status', 'Date', 'IPR_annotation', 'IPR_description', 'GO_annotation',
               'Pathway_annotation']
    if os.path.getsize(outfile_Interpro):
        with open(outfile_Interpro) as handle:
            Interpro_df = pandas.read_csv(handle, sep='\t', index_col=None, names=headers, converters={'ID': strip})

        Interpro_df.columns = Interpro_df.columns.map(lambda x: str(x) + '_Interpro')

        Interpro_reduced_df = Interpro_df.groupby('ID_Interpro', as_index=False)['IPR_annotation_Interpro', 'IPR_description_Interpro'].agg(lambda x: list(x))
        results_df[['Domain_type', 'IPR_annotation_Interpro', 'IPR_description_Interpro']] = results_df.apply(lambda row:
                                                pandas.Series(Interpro_arrange(row, Domains, Interpro_reduced_df, 'IPR')), axis=1)
        ni_reduced_df = Interpro_df.groupby('ID_Interpro', as_index=False)['Signature Accession_Interpro', 'Signature Description_Interpro'].agg(lambda x: list(x))
        results_df[['Domain_type_ni', 'ni_annotation_Interpro', 'ni_description_Interpro']] = results_df.apply(lambda row:
                                                pandas.Series(Interpro_arrange(row, Domains_ni, ni_reduced_df, 'NI')), axis=1)
        go_reduced_df = Interpro_df.groupby('ID_Interpro', as_index=True)['GO_annotation_Interpro'].agg(lambda x: list(x))
        go_reduced_df = go_reduced_df.to_frame().reset_index()
        go_reduced_df['GO_annotation_Interpro'] = go_reduced_df.apply(lambda row: split_GO(row), axis=1)
        results_df[['GO_annotation_Interpro','GO_terms_Interpro']] = results_df.apply(lambda row:
                                            pandas.Series(Interpro_arrange(row, GO_dict, go_reduced_df, 'GO')), axis=1)

    else:
        cols_ipr_go = ['Domain_type', 'IPR_annotation_Interpro',
                    'IPR_description_Interpro', 'Domain_type_ni',
                    'ni_annotation_Interpro', 'ni_description_Interpro',
                    'GO_annotation_Interpro', 'GO_terms_Interpro']
        results_df = pandas.concat([results_df, pandas.DataFrame(columns=cols_ipr_go)], axis=1, sort=False)


    results_df[['Predicted_SCL', 'Notes_SCL']] = results_df.apply(lambda row: pandas.Series(Predict(row)), axis=1)

    results_df = results_df.drop(labels=['# ID_SigP5', 'SEQENCE_ID_Phobius', 'FASTA-ID_PrediSi', 'ID_TatP', 'ID_LipoP', 'ID_TMHMM', 'ID_PCB'], axis=1)
    results_df = results_df.rename(columns={"name_SigP4": "Protein_ID"})

    if args.run == 'all':
        results_df['Protein_ID'] = old_ids
    elif args.run == 'pred':
        with open(output_path + new_foldername + '/ID_' + new_foldername + '.id', 'r') as handle:
            old_ids = handle.read().splitlines()
        results_df['Protein_ID'] = old_ids

    with open(outfile_results, 'w') as handle:
        results_df.to_csv(handle, sep='\t')



t_end = datetime.datetime.now()
t_diff = t_end - t_start
logging.info("Run finished at " + str(t_end) + ", in " + str(t_diff) + ". You can find the final results in: " + outfile_results)
if args.run == 'data' or args.run == 'all':
    shutil.move('GP4_'+log_name+'.log', output_path + new_foldername + '/')
elif args.run == 'pred':
    shutil.move(script_path +'/' + 'GP4_' + log_name + '.log', output_path + new_foldername + '/' +'GP4_' + log_name + '_2' + '.log')
