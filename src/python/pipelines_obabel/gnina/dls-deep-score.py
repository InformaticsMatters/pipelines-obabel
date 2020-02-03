#!/usr/bin/env python

# Copyright 2020 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# Create dir containing ligands.sdf and protein.pdb
# Enter docker container like this:
#   docker run -it --rm --gpus all -v $PWD:/root/train/fragalysis_test_files/work:Z informaticsmatters/deep-app-ubuntu-1604:latest bash
#
# Now inside the container run like this:
#   rm -rf work/* && python3 -m  pipelines_obabel.gnina.dls-deep-score -i ligands.sdf -r protein.pdb -w work
#
# If testing with no GPU you can use the --mock option to generate random scores
#
# Start container for testing like this:
#    docker run -it --rm -v $PWD:$PWD:Z -w $PWD --entrypoint bash pipelines-obabel-dls-deep
# Inside container test like this:
#   mkdir /tmp/work
#   rm -rf /tmp/work/* && python3 -m pipelines_obabel.gnina.dls-deep-score -i ../../data/ligands.sdf.gz -r ../../data/nudt7/receptor.pdb.gz -w /tmp/work --mock
#

import argparse, os, re
import random
from openbabel import pybel
from pipelines_obabel import obabel_utils
from pipelines_utils import utils, parameter_utils

types_file_name = 'inputs.types'
types_file_name = 'inputs.types'
predict_file_name = 'predictions.txt'
work_dir = '.'
inputs_protein = []
inputs_ligands = []

def write_inputs(protein_file, ligands_file):
    global types_file_name
    global work_dir
    global inputs_protein
    global inputs_ligands

    ligands_path = "{0}{1}ligands".format(work_dir, os.path.sep)
    utils.log("Writing ligands to", ligands_path)
    os.mkdir(ligands_path)
    cmd1 = "gninatyper {0} {1}{2}ligands{2}ligand".format(ligands_file, work_dir, os.path.sep)
    utils.log('CMD:', cmd1)
    os.system(cmd1)
    ligand_gninatypes = os.listdir("{0}{1}ligands".format(work_dir, os.path.sep))

    proteins_path = "{0}{1}proteins".format(work_dir, os.path.sep)
    utils.log("Writing proteins to", proteins_path)
    os.mkdir(proteins_path)
    cmd2 = "gninatyper {0} {1}{2}proteins{2}protein".format(protein_file, work_dir, os.path.sep)
    utils.log('CMD:', cmd2)
    os.system(cmd2)
    protein_gninatypes = os.listdir("{0}{1}proteins".format(work_dir, os.path.sep))

    types_path = "{0}{1}{2}".format(work_dir, os.path.sep, types_file_name)
    utils.log("Writing types to", types_path)
    with open(types_path, 'w') as types_file:
        for protein in protein_gninatypes:
            inputs_protein.append(protein)
            for ligand in ligand_gninatypes:
                utils.log("Handling", protein, ligand)
                inputs_ligands.append(ligand)
                line = "0 {0}{3}proteins{3}{1} {0}{3}ligands{3}{2}\n".format(work_dir, protein, ligand, os.path.sep)
                types_file.write(line)

def generate_predictions_filename(work_dir, predict_file_name):
    return "{0}{1}{2}".format(work_dir, os.path.sep, predict_file_name)

def run_predictions():
    global types_file_name
    global predict_file_name
    global work_dir
    # python3 scripts/predict.py -m resources/dense.prototxt -w resources/weights.caffemodel -i work_0/test_set.types >> work_0/caffe_output/predictions.txt
    cmd1 = "python3 /root/train/fragalysis_test_files/scripts/predict.py -m /root/train/fragalysis_test_files/resources/dense.prototxt" +\
           " -w /root/train/fragalysis_test_files/resources/weights.caffemodel" +\
            " -i {0}/{1} -o {0}/{2}".format(work_dir, types_file_name, predict_file_name)
    utils.log("CMD:", cmd1)
    os.system(cmd1)

def mock_predictions():
    global work_dir
    global predict_file_name
    global inputs_protein
    global inputs_ligands
    utils.log("WARNING: generating mock results instead of running on GPU")
    outfile = generate_predictions_filename(work_dir, predict_file_name)
    with open(outfile, 'w') as predictions:
        for protein in inputs_protein:
            for ligand in inputs_ligands:
                score = random.random()
                line = "{0} | 0 {1}{4}proteins{4}{2} {1}{4}ligands{4}{3}\n".format(score, work_dir, protein, ligand, os.path.sep)
                # utils.log("Writing", line)
                predictions.write(line)


def read_predictions():
    global predict_file_name
    global work_dir
    scores = {}
    with open("{0}{1}{2}".format(work_dir, os.path.sep, predict_file_name), 'r') as input:
        for line in input:
            #utils.log(line)
            tokens = line.split()
            if len(tokens) == 5 and tokens[1] == '|':
                # utils.log(len(tokens), tokens[0], tokens[3], tokens[4])
                record_no = match_ligand(tokens[4])
                if record_no is not None:
                    # utils.log(record_no, tokens[0])
                    scores[record_no] = tokens[0]
    utils.log("Found", len(scores), "scores")
    return scores

def patch_scores_sdf(sdf_in, outfile, scores):

    global work_dir

    counter = 0
    sdf_path = "{0}{1}{2}.sdf".format(work_dir, os.path.sep, outfile)
    tsv_path = "{0}{1}{2}.tsv".format(work_dir, os.path.sep, outfile)
    utils.log("Writing results to {0} and {1}".format(sdf_path, tsv_path))
    with open(tsv_path, 'w') as tsv_file:
        sdf_file = pybel.Outputfile("sdf", sdf_path)
        for mol in pybel.readfile("sdf", sdf_in):
            if counter in scores:
                score = scores[counter]
                # utils.log("Score for record {0} is {1}".format(counter, score))

                mol.data['dls_deep_score'] = score
                if 'SCORE' in mol.data:
                    rdock_score = mol.data['SCORE']
                else:
                    rdock_score = ''

                if 'SCORE.norm' in mol.data:
                    rdock_nscore = mol.data['SCORE.norm']
                else:
                    rdock_nscore = ''

                sdf_file.write(mol)
                tsv_file.write("{0}\t{1}\t{2}\t{3}\n".format(counter, rdock_score, rdock_nscore, score))

            else:
                utils.log("No score found for record", counter)
            counter += 1
        sdf_file.close()

def patch_scores_json(sdf_in, outfile, scores):

    global work_dir

    counter = 0
    json_path = "{0}{1}{2}".format(work_dir, os.path.sep, outfile)
    output = utils.open_output(json_path, 'data', True)
    writer = obabel_utils.json_writer(output)

    tsv_path = "{0}{1}{2}.tsv".format(work_dir, os.path.sep, outfile)
    utils.log("Writing results to {0} and {1}".format(output.name, tsv_path))
    with open(tsv_path, 'w') as tsv_file:
        for mol in pybel.readfile("sdf", sdf_in):
            if counter in scores:
                score = scores[counter]
                # utils.log("Score for record {0} is {1}".format(counter, score))

                mol.data['dls_deep_score'] = score
                if 'SCORE' in mol.data:
                    rdock_score = mol.data['SCORE']
                else:
                    rdock_score = ''

                if 'SCORE.norm' in mol.data:
                    rdock_nscore = mol.data['SCORE.norm']
                else:
                    rdock_nscore = ''

                writer.write(mol)
                tsv_file.write("{0}\t{1}\t{2}\t{3}\n".format(counter, rdock_score, rdock_nscore, score))

            else:
                utils.log("No score found for record", counter)
            counter += 1
            writer.flush()
        writer.close()

# work/ligands/ligand_9.gninatypes
ligand_patt = re.compile(r'.*ligands/ligand_(\d+)\.gninatypes$')

def match_ligand(s):
    global ligand_patt
    m = ligand_patt.match(s)
    if m:
        i = m.group(1)
        return int(i)
    else:
        return None

def write_json_as_sdf(jsonfile, sdfile):
    input, suppl = obabel_utils.default_open_input(jsonfile, 'json')
    with pybel.Outputfile('sdf', sdfile) as output:
        for mol in suppl:
            output.write(mol)


def main():

    global work_dir


    parser = argparse.ArgumentParser(description='DLS Deep - pose scoring')
    parameter_utils.add_default_input_args(parser)
    parser.add_argument('--no-gzip', action='store_true', help='Do not compress the output (STDOUT is never compressed')
    parser.add_argument('-r', '--receptor', help="Receptor file for scoring (PDB or Mol2 format)")
    parser.add_argument('-o', '--outfile', default='scored_ligands', help="Base file name for results")
    parser.add_argument('-of', '--outformat', choices=['sdf', 'json'], default='sdf', help="Output format. Defaults to 'sdf'.")
    parser.add_argument('-w', '--work-dir', default=".", help="Working directory")
    parser.add_argument('--mock', action='store_true', help='Generate mock scores rather than run on GPU')
    parser.add_argument('--thin', action='store_true', help='Thin output mode')

    args = parser.parse_args()
    utils.log("DLS deep args: ", args)

    work_dir = args.work_dir

    informat = args.informat
    protein = args.receptor
    ligands = args.input
    outfile = args.outfile

    if informat == 'json' or ligands.lower().endswith('.data') or ligands.lower().endswith('.data.gz'):
        # we need to write to SDF
        utils.log("Converting ligands from JSON to SDF")
        ligands_sdf = "{0}{1}ligands.sdf".format(work_dir, os.path.sep)
        write_json_as_sdf(ligands, ligands_sdf)

    elif informat == 'sdf' or ligands.lower().endswith('.sdf') or ligands.lower().endswith('.sdf.gz'):
        ligands_sdf = ligands


    else:
        raise ValueError("Unexpected input format for ligands")


    # # Open the output file
    # s_now = datetime.datetime.utcnow().strftime("%d-%b-%Y %H:%M:%S UTC")
    # source = 'pipelines/gnina/dls-deep-score.py'
    # output, WRITER, output_base = \
    #     rdkit_utils.default_open_output(args.output, "dls-deep-score", args.outformat,
    #                                     compress=not args.no_gzip,
    #                                     thinOutput=args.thin,
    #                                     valueClassMappings={'dls-deep-score': 'java.lang.Float'},
    #                                     datasetMetaProps={'created': s_now,
    #                                                       'source': source,
    #                                                       'description': 'DLS Deep - pose scoring'}
    #                                     )
    #
    # PDB_PATH = args.pdb_file
    #    # Close the file
    # WRITER.close()

    write_inputs(protein, ligands_sdf)
    if args.mock:
        mock_predictions()
    else:
        run_predictions()
    scores = read_predictions()

    if args.outformat == 'sdf':
        patch_scores_sdf(ligands_sdf, outfile, scores)
    elif args.outformat == 'json':
        patch_scores_json(ligands_sdf, outfile, scores)

    if args.outformat == 'sdf':
        if not args.no_gzip:
            os.system("gzip {0}{1}{2}.sdf".format(work_dir, os.path.sep, outfile))

    # if args.meta:
    #     utils.write_metrics(output_base, {'__InputCount__': COUNTER, '__OutputCount__': SUCCESS, 'PLI': COUNTER})


if __name__ == "__main__":
    main()
