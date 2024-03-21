import requests
import codecs
import json
import urllib
from .AnnotatedSpectrum import *
import sys


class XiAnnotatorSuper:
    def __init__(self):
        self.last_request = None
        self.last_response = None
        self.annotated_spectrum = AnnotatedSpectrum()

    @staticmethod
    def create_json_annotation_request(
            peak_list,
            peptide,
            precursor_charge,
            precursor_mz,
            fragment_types=('peptide', 'b', 'y'),
            fragment_tolerance_ppm=0.0,
            fragment_tolerance_abs=0.0,
            cross_linker="BS3",
            custom_settings=False,
            as_dict=False,
            modifications=None
    ):
        """
        Runs the Xi annotator to get the fragment annotation.
        :param peak_list: list of lists [[mz, int], ...]
        :param peptide: Peptide Object
        :param precursor_charge:
        :param precursor_mz:
        :param fragment_types: ['peptide', 'a', 'b', 'c', 'x', 'y', 'z']
        :param fragment_tolerance_ppm: (float) fragment tolerance in ppm
        :param fragment_tolerance_abs: (float) fragment tolerance in Dalton
        :param cross_linker: either (str) cross-linker used or (float) cross-linker mod mass
        :param custom_settings: (list) custom_settings
        :param as_dict: (bool) returns request as dictionary instead of json
        :param modifications: (dict) define additional modifications as name:mass dict
        :return: annotation request
        """

        if modifications is None:
            modifications = {}
        if type(cross_linker) is float:
            cross_linker_mod_mass = cross_linker
        else:
            cross_linker_mod_masses = {
                'BS3': 138.06807961,
                'DSSO': 158.0037648,
                "UCL2": 250.012224,
                "UCL3": 218.040152
            }
            try:
                cross_linker_mod_mass = cross_linker_mod_masses[cross_linker]
            except KeyError:
                raise ValueError('unknown crosslinker: %s' % cross_linker)

        mod_mass_dict = {
            'bs3': 27.984,  # ???
            'bs3loop': 138.06807,
            'bs3nh': 155.094619105,
            'bs3nh2': 155.094619105,
            'bs3oh': 156.0786347,
            'cm': 57.021464,
            'ox': 15.994915,
            'dsso': 158.0037648,
            'dssonh2': 175.030313905,
            'dssooh': 176.0143295
        }
        mod_mass_dict.update(modifications)

        if fragment_tolerance_ppm != 0 and fragment_tolerance_abs != 0:
            raise ValueError('Give either fragment_tolerance_ppm or fragment_tolerance_abs')
        if fragment_tolerance_ppm != 0:
            frag_tolerance_val = str(fragment_tolerance_ppm)
            frag_tolerance_unit = 'ppm'
        else:
            frag_tolerance_val = str(fragment_tolerance_abs)
            frag_tolerance_unit = 'Da'

        all_mods = []

        # peptides block
        peptides = []

        pep1 = [{'aminoAcid': char, 'Modification': ''} for char in peptide.pep_seq1 if char.isupper()]
        offset = 1
        for match in re.finditer('[^A-Z]+', peptide.pep_seq1):
            modification = match.group()
            pep1[match.start() - offset]['Modification'] = modification
            offset += len(modification)
            # add to all mods
            if modification not in all_mods:
                all_mods.append(modification)

        peptides.append(pep1)

        if not peptide.is_linear:
            pep2 = [{'aminoAcid': char, 'Modification': ''} for char in peptide.pep_seq2 if char.isupper()]
            offset = 1
            for match in re.finditer('[^A-Z]+', peptide.pep_seq2):
                modification = match.group()
                pep2[match.start() - offset]['Modification'] = modification
                offset += len(modification)
                # add to all mods
                if modification not in all_mods:
                    all_mods.append(modification)

            peptides.append(pep2)

        peptides_block = [{"sequence": pep} for pep in peptides]

        # link sites block
        if peptide.is_linear:
            link_sites_block = []
        else:
            link_sites_block = [
                {'id': 0, 'peptideId': 0, 'linkSite': peptide.link_pos1 - 1},
                {'id': 0, 'peptideId': 1, 'linkSite': peptide.link_pos2 - 1}
            ]

        # peak list
        peak_block = [{"mz": float(x[0]), "intensity": float(x[1])} for x in peak_list]

        # annotation block
        annotation_modifications = [{"aminoAcids": ["*"], "id": mod, "mass": mod_mass_dict[mod]} for mod in all_mods]

        ion_types = [{'type': ion.title() + 'Ion'} for ion in fragment_types]

        annotation_block = {
            "fragmentTolerance": {"tolerance": frag_tolerance_val, "unit": frag_tolerance_unit},
            "modifications": annotation_modifications,
            "ions": ion_types,
            "crosslinker": {"modMass": cross_linker_mod_mass},
            "precursorCharge": precursor_charge,
            "precursorMZ": precursor_mz,
            "custom": custom_settings
        }

        json_dict = {
            "Peptides": peptides_block,
            "LinkSite": link_sites_block,
            "peaks": peak_block,
            "annotation": annotation_block
        }

        if as_dict:
            return json_dict

        return json.dumps(json_dict)

    def get_annotated_spectrum(self):
        return self.annotated_spectrum

    def get_json_response(self):
        return self.last_response

    def get_json_request(self):
        return self.last_request


class XiAnnotatorLocal(XiAnnotatorSuper):
    def __init__(self, java_home_dir='/usr/lib/jvm/java-8-openjdk-amd64/', jar_path='', jvmargs=''):
        XiAnnotatorSuper.__init__(self)

        self.xiAnnotatorVersion = '1.4.28'

        import os
        os.environ['JAVA_HOME'] = java_home_dir

        if jar_path == '':
            if getattr(sys, 'frozen', False):
                jar_root = os.path.dirname(sys.executable)
            else:
                jar_root = os.path.dirname(__file__)

            jar_path = os.path.join(
                jar_root,
                'jar',
                'xiAnnotator-{0}-jar-with-dependencies.jar'.format(self.xiAnnotatorVersion)
            )

        if 'CLASSPATH' in os.environ.keys():

            jars = os.environ['CLASSPATH'].split(os.pathsep)
            jars = [j for j in jars if 'xiAnnotator' not in j]
            jars.append(jar_path)
            os.environ['CLASSPATH'] = os.pathsep.join(jars)
        else:
            os.environ['CLASSPATH'] = jar_path

        import jpype
        if not jpype.isJVMStarted():
            jpype.startJVM(jvmargs, classpath=jar_path)
        self.JString = jpype.java.lang.String
        self.annotator = jpype.JPackage('org').rappsilber.xiAnnotator()

    def request_annotation_json(self, json_request):
        """
        :param json_request: xiAnnotator json request
        :return: annotation response JSON
        """

        self.last_request = json_request

        response = self.annotator.getFullAnnotation(self.JString(json_request)).getEntity()
        json_response = json.loads(str(response))

        self.last_response = json_response
        self.annotated_spectrum.load_json(json_response)


class XiAnnotatorWeb(XiAnnotatorSuper):
    def __init__(self, url):
        self.base_url = url
        XiAnnotatorSuper.__init__(self)

    def get_known_modifications(self):
        """
        Returns a list of known modifications
        in the form of a dict: aminoAcids: [char]
                              id: str
                              mass: float
        """
        url = self.base_url + '/annotate/knownModifications'
        response = requests.get(url)
        if response.status_code == 200:
            return json.loads(response.content)['modifications']
        else:
            raise Exception(response.status_code)

    @staticmethod
    def request_annotation_xidb(search_id, psm_id, pep_seq1, pep_seq2, link_pos1, link_pos2):
        """
        Runs the Xi annotator to get the fragment annotation.
        """
        if search_id > 10000:
            base_url = 'http://xi3.bio.ed.ac.uk/xiAnnotator/annotate/'
        else:
            base_url = 'http://129.215.14.63/xiAnnotator/annotate/'

        url = base_url + '%s/85160-94827-76653-69142/%s/?peptide=%s&peptide=%s&link=%s&link=%s' % (
            int(search_id), int(psm_id), pep_seq1, pep_seq2, int(link_pos1), int(link_pos2))

        reader = codecs.getreader("utf-8")
        data = json.load(reader(urllib.request.urlopen(url)))
        return data

    def request_annotation_json(self, json_request):
        """

        :param json_request: xiAnnotator json request
        :return: annotation response JSON
        """
        headers = {'Content-type': 'application/json', 'Accept': 'application/json'}
        url = self.base_url + 'annotate/FULL'
        r = requests.post(url, data=json_request, headers=headers)
        self.last_request = json_request
        response_json = r.json()
        self.last_response = response_json
        self.annotated_spectrum.load_json(response_json)
