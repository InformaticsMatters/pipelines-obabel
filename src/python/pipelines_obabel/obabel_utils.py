#!/usr/bin/env python

# Copyright 2019 Informatics Matters Ltd.
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


import sys, gzip, json, uuid
from openbabel import pybel

from pipelines_utils import utils
from pipelines_utils.StreamJsonListLoader import StreamJsonListLoader


def default_open_input(inputDef, inputFormat):
    if not inputDef and not inputFormat:
        raise ValueError('Must specify either an input file name or an input format (or both)')
    elif inputFormat == 'sdf' or (inputDef and (inputDef.lower().endswith('.sdf') or inputDef.lower().endswith('.sdf.gz'))):
        input, suppl = default_open_input_sdf(inputDef)
    elif inputFormat == 'json' or (inputDef and (inputDef.lower().endswith('.data') or inputDef.lower().endswith('.data.gz'))):
        input, suppl = default_open_input_json(inputDef)
    # elif inputFormat == 'smiles':
    #     input, suppl = default_open_input_typed_smiles(inputDef)
    else:
        raise ValueError('Unsupported input format')

    return input, suppl


def default_open_input_sdf(filename):
    """Open the input as a SD file (possibly gzipped if ending with .gz) according to RDKit's ForwardSDMolSupplier

    :param filename: The name of the file.
    """

    suppl = pybel.readfile("sdf", filename)
    return None, suppl

def default_open_input_json(filename, lazy=True):
    """Open the given input as JSON array of Squonk MoleculeObjects.

    :param inputDef: The name of the input file, or None if to use STDIN.
                     If filename ends with .gz will be gunzipped
    :param lazy: Use lazy loading of the JSON. If True will allow handling of
                 large datasets without being loaded into memory,
                 but may be less robust and will be slower.
    """

    if filename.lower().endswith('.gz'):
        input = gzip.open(filename, 'rt')
    else:
        input = open(filename, 'r')

    if lazy:
        suppl = generate_mols_from_json(StreamJsonListLoader(input))
    else:
        suppl = generate_mols_from_json(json.load(input))
    return input, suppl

def generate_mols_from_json(input):
    """Create a supplier of OBMol Mol objects from the json

    :param input: file like object containing the json representation of the molecules
    """
    j=0
    for item in input:
        j+=1
        mol = create_mol_from_props(item)
        if not mol:
            # TODO - get a count of the errors and report
            utils.log("Failed to create molecule - skipping. Data was ", item)
            continue
        yield mol


def create_mol_from_props(molobj):
    """Function to get the RDKit mol from MoleculeObject JSON

    :param molobj: Python dictionary containing the molecule's properties
    """

    if "source" not in molobj or "format" not in molobj:
        return None

    molstr = str(molobj["source"])
    # Get the format and use this as a starting point to work out
    molformat = molobj["format"]
    # Now parse it with obabel
    mol = parse_mol_simple(molformat, molstr)
    if mol:
        if "values" in molobj:
            values = molobj["values"]
            for key in values:
                mol.data[str(key)] = values[key]
        uuid = str(molobj["uuid"])
        if uuid:
            mol.data["uuid"] = uuid
            mol.data["_Name"] = uuid
    return mol

def parse_mol_simple(molformat, molstr):
    if molformat == "smiles":
        format = "smi"
    else:
        format = molformat

    return pybel.readstring(format, molstr)

def default_open_output_json(outputDef, outputBase, thinOutput,
                             compress, valueClassMappings, datasetMetaProps,
                             fieldMetaProps):

    # this writes the metadata that Squonk needs
    utils.write_squonk_datasetmetadata(outputBase, thinOutput, valueClassMappings,
                                       datasetMetaProps, fieldMetaProps)

    output = utils.open_output(outputDef, 'data', compress)

    if thinOutput:
        # writer = ThinJsonWriter(output)
        raise("Thin output not yet supported")
    else:
        writer = ThickJsonWriter(output)
    return output,writer

def json_writer(output):
    return ThickJsonWriter(output)

class ThickJsonWriter:

    def __init__(self, file):
        self.file = file
        self.file.write('[')
        self.count = 0

    def write(self, mol, props=None, includeStereo=True, confId=-1,
              kekulize=True, forceV3000=False, format='mol'):
        d = {}
        if format == 'mol':
            # d['source'] = Chem.MolToMolBlock(mol, includeStereo=includeStereo, confId=confId, kekulize=kekulize, forceV3000=forceV3000)
            s = mol.write(format='mol')
            # utils.log("Writing mol as", s)
            d['source'] = s
            d['format'] = 'mol'
        elif format == 'smiles':
            # if kekulize:
            #     Chem.Kekulize(mol)
            d['source'] = mol.write(format='smi')
            d['format'] = 'smiles'
        else:
            raise ValueError("Unexpected format: " + format)
        allProps = {}
        for name in mol.data:
            if name not in ['MOL Chiral Flag', 'OpenBabel Symmetry Classes']:
                allProps[name] = mol.data[name]

        if props:
            allProps.update(props)

        # print("Props: " + ",".join(allProps))

        if 'uuid' in allProps:
            d['uuid'] = allProps['uuid']
            del allProps['uuid']
        else:
            d['uuid'] = str(uuid.uuid4())
        if allProps:
            d['values'] = allProps
        #utils.log("Mol:",d)
        json_str = json.dumps(d)
        if self.count > 0:
            self.file.write(',')
        self.file.write(json_str)
        self.count += 1

    def close(self):
        self.file.write(']')
        self.file.close()

    def flush(self):
        self.file.flush()