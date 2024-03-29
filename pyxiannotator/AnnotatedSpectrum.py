import re
import math
from memoized_property import memoized_property
import numpy as np


class AnnotatedSpectrum:
    def __init__(self):
        self.peaks = []
        self.clusters = []
        self.deisotoped_sum_peaks = []
        self.deisotoped_max_peaks = []
        self.fragments = []
        self.annotation_json = None
        self.precursor = {}
        self.isLinear = None
        self.peptide = None

        self.base_peak = None
        self.deisotoped_sum_base_peak = None
        self.deisotoped_max_base_peak = None
        self.intensity_sorted_peaks = None
        self.int_sorted_deisotoped_sum_peaks = None
        self.int_sorted_deisotoped_max_peaks = None

        self.deisotoped_sum_ranks = None
        self.deisotoped_max_ranks = None
        self.ranks = None

    def _deisotope_peaks_(self, deisotope_func):
        """
        Deisotope the peaks.

        :param deisotope_func: (str) function to use for deisotoped intensity.
            Valid values: - 'sum': summed up cluster intensity
                          - 'max': maximum peak intensity in the cluster
        """
        non_cluster_peaks = [p for p in self.peaks if len(p.cluster_ids) == 0]

        cluster_peaks = [
            Peak(cluster.get_first_peak().id, cluster.get_first_peak().mz, cluster.get_intensity(deisotope_func),
                 self, [cluster.id])
            for cluster in self.clusters]

        deisotoped_peaks = non_cluster_peaks + cluster_peaks

        return sorted(deisotoped_peaks, key=lambda k: k.mz)

    def load_json(self, annotation_json):

        def create_fragments(anno_json):
            found_fragments = []

            for fragment in anno_json['fragments']:
                # each cluster is one charge state
                for fragment_cluster_info in fragment['clusterInfo']:
                    # assert fragment_cluster_info['errorUnit'] == 'ppm'

                    fragment_cluster_id = fragment_cluster_info['Clusterid']
                    # fragment_cluster = anno_json['clusters'][fragment_cluster_id]

                    # fragment_cluster_peaks = [
                    #     Peak(peak['mz'], peak['intensity']) for peak in anno_json['peaks'] if
                    #     fragment_cluster_id in peak['clusterIds']]

                    miss_monoiso = bool(fragment_cluster_info['matchedMissingMonoIsotopic'])

                    if fragment['class'] == 'lossy':
                        lossy = True
                    elif fragment['class'] == 'non-lossy':
                        lossy = False
                    else:
                        raise Exception('Unknown fragment class: %s' % fragment['class'])

                    error = {'value': fragment_cluster_info['error'],
                             'unit': fragment_cluster_info['errorUnit']}

                    found_fragments.append(
                        Fragment(
                            fragment['name'],
                            fragment['peptideId'],
                            self.clusters[fragment_cluster_id],
                            # IsotopeCluster(fragment_cluster_peaks, fragment_cluster['charge']),
                            lossy,
                            miss_monoiso,
                            fragment['sequence'],
                            fragment_cluster_info['calcMZ'],
                            error,
                            self.clusters[fragment_cluster_id].get_first_peak(),
                            self,
                            fragment['type'],
                            fragment_cluster_info['matchedCharge']
                        )
                    )
            return found_fragments

        self.peaks = [Peak(i, p['mz'], p['intensity'], self, p['clusterIds'])
                      for i, p in enumerate(annotation_json['peaks'])]

        self.clusters = []
        cluster_index = 0
        for i, cluster in enumerate(annotation_json['clusters']):
            cluster_peaks = [p for p in self.peaks if cluster_index in p.cluster_ids]
            self.clusters.append(IsotopeCluster(i, cluster_peaks, cluster['charge']))
            cluster_index += 1

        self.deisotoped_sum_peaks = self._deisotope_peaks_(deisotope_func='sum')
        self.deisotoped_max_peaks = self._deisotope_peaks_(deisotope_func='max')

        # get base peak (deisotoped and non deisotoped)
        self.deisotoped_sum_base_peak = sorted(self.deisotoped_sum_peaks,
                                               key=lambda k: k.intensity)[-1]
        self.deisotoped_max_base_peak = sorted(self.deisotoped_max_peaks,
                                               key=lambda k: k.intensity)[-1]
        self.base_peak = sorted(self.peaks, key=lambda k: k.intensity)[-1]

        # intensity sorted peaks
        deisotoped_sum_peaks = self.get_peaks(deisotoped=True, as_list=True, deisotope_func='sum')
        self.int_sorted_deisotoped_sum_peaks = sorted(deisotoped_sum_peaks, key=lambda p: p[1])
        deisotoped_max_peaks = self.get_peaks(deisotoped=True, as_list=True, deisotope_func='max')
        self.int_sorted_deisotoped_max_peaks = sorted(deisotoped_max_peaks, key=lambda p: p[1])
        peaks = self.get_peaks(deisotoped=False, as_list=True)
        self.intensity_sorted_peaks = sorted(peaks, key=lambda p: p[1])

        # ranks
        def calc_ranks(peak_list):
            """
            Calculate ranks of peaks (same intensity = same rank)
            :param peak_list: list of Peak
            :return: (ndarray) ranks
            """
            intensities = np.array([p.intensity for p in peak_list])
            sorted_indices = np.argsort(intensities)
            sorted_intensities = intensities[sorted_indices]
            unique_intensities, indices, ranks = np.unique(sorted_intensities, return_index=True,
                                                           return_inverse=True)
            # ranks are the wrong way round - change this and convert to 1-based
            ranks = ranks.max(initial=0) - ranks + 1
            return ranks

        self.deisotoped_sum_ranks = calc_ranks(self.deisotoped_sum_peaks)
        self.deisotoped_max_ranks = calc_ranks(self.deisotoped_max_peaks)
        self.ranks = calc_ranks(self.peaks)

        self.fragments = create_fragments(annotation_json)

        self.annotation_json = annotation_json

        self.precursor = {
            'mz': annotation_json['annotation']['precursorMZ'],
            'intensity': annotation_json['annotation']['precursorIntensity'],
            'charge': annotation_json['annotation']['precursorCharge'],
        }

        pep_seqs = []
        for peptide_json in annotation_json['Peptides']:
            pep_seq = "".join(
                [aa['aminoAcid'] + aa['Modification'] for aa in peptide_json['sequence']])
            pep_seqs.append(pep_seq)

        if len(pep_seqs) == 1:
            pep_seqs.append('')
            self.isLinear = True

        try:
            # link_pos in annotator is 0 based
            link_pos1 = annotation_json['LinkSite'][0]['linkSite'] + 1
            link_pos2 = annotation_json['LinkSite'][1]['linkSite'] + 1
        except IndexError:
            link_pos1 = -1
            link_pos2 = -1

        self.peptide = MzSpecies(
            pep_seqs[0],
            pep_seqs[1],
            annotation_json['annotation']['precursorCharge'],
            None,  # annotation_json['annotation']['precursorMZ'],
            # ToDo: precursorMZ = matchMZ? should maybe use calcMZ?
            link_pos1,
            link_pos2
        )

    def get_fragments(self, lossy=False, as_dict=False):
        """
        get the fragments of the spectrum. Either as list of Fragment objects or list of dictionaries
        :param lossy: include neutral loss fragments
        :param as_dict: returns list of dictionaries
        :return: matched fragments
        """
        if lossy:
            fragments = self.fragments
        else:
            fragments = [f for f in self.fragments if not f.lossy]

        if as_dict:
            return [f.as_dict() for f in fragments]
        else:
            return fragments

    def get_peaks(self, as_list=False, deisotoped=False, deisotope_func='sum'):
        """
        Get the peaks of the spectrum.

        Either as list of Peak objects or list of [mz, int] with the option of deisotoping according
         to found clusters.
        :param as_list: return as list of [mz, int]
        :param deisotoped: deisotoping on/off
        :param deisotope_func: (str) function to use for deisotoped intensity.
            Valid values: - 'sum': summed up cluster intensity
                          - 'max': maximum peak intensity in the cluster
        :return: spectrum peaks
        """
        if deisotoped:
            if deisotope_func == 'sum':
                peaks = self.deisotoped_sum_peaks
            elif deisotope_func == 'max':
                peaks = self.deisotoped_max_peaks
            else:
                raise ValueError(f'Unknown deisotope_func: {deisotope_func}')
        else:
            peaks = self.peaks

        if as_list:
            return [p.as_list() for p in peaks]

        return peaks

    def get_base_peak(self, as_list=False, deisotoped=False):
        """
        get the most intense peak of the spectrum
        :param as_list: return as list of [mz, int]
        :param deisotoped: deisotoping on/off
        :return: base peak
        """
        if deisotoped:
            # Fixme
            raise NotImplementedError('deisotoped base peak not implemented')
            # base_peak = self.deisotoped_base_peak
        else:
            base_peak = self.base_peak

        if as_list:
            return base_peak.as_list()
        else:
            return base_peak

    def get_clusters(self):
        return self.clusters

    def get_cluster_by_id(self, cluster_id):
        return self.clusters[cluster_id]

    def match_peak(self, needle_mz, tolerance, deisotoped=False):
        """
        Find a peak in spectrum with a certain tolerance.

        (default fragment tolerance from annotation)
        :param needle_mz: needle m/z value
        :param tolerance: relative error tolerance for matching
        :param deisotoped: deisotoping on/off
        :return: peak, error
        """
        peaks = self.get_peaks(deisotoped=deisotoped)

        errors = [math.fabs(haystack_peak.match_error(needle_mz)) for haystack_peak in peaks]

        abs_match_error = min(errors)

        if abs_match_error > tolerance:
            return None, None

        match_peak = peaks[errors.index(abs_match_error)]

        match_err = (needle_mz - match_peak.mz) / match_peak.mz

        return match_peak, match_err

    def match_unfragmented_precursor_peak(self, tolerance, deisotoped=False):
        """
        :param tolerance: error tolerance for matching (default fragment tolerance from annotation)
        :param deisotoped: deisotoping on/off
        :return: precursor peak, error
        """
        precursor_mz = self.precursor['mz']

        return self.match_peak(precursor_mz, tolerance, deisotoped)

    def get_unfragmented_precursor_fragments(self):
        """
        Return the unfragmented precursor (matched precursor fragments).

        :return: unfragmented precursor fragment
        """
        precursor_fragments = [
            f for f in self.get_fragments() if f.ion_type == 'Precursor'
            # and f.get_charge() == self.precursor['charge']
        ]
        if len(precursor_fragments) == 0:
            return None

        return precursor_fragments

    def get_unfragmented_precursor_intensity(self, deisotoped=False):
        """Return the unfragmented precursor intensity."""
        precursor_frags = self.get_unfragmented_precursor_fragments()
        if precursor_frags is None:
            return 0
        if deisotoped:
            cluster = [f.cluster for f in precursor_frags]
            peaks = [c.peaks for c in cluster]
            peaks = [item for sublist in peaks for item in sublist]
        else:
            peaks = [f.peak for f in precursor_frags]
        unique_peak_ids = set([p.id for p in peaks])
        summed_intensity = sum([self.peaks[j].intensity for j in unique_peak_ids])
        return summed_intensity

    def get_peak_rank(self, peak, deisotoped=False, deisotope_func='sum', as_list=False):
        if deisotoped:
            if deisotope_func == 'sum':
                peaks = self.int_sorted_deisotoped_sum_peaks
                ranks = self.deisotoped_sum_ranks
            elif deisotope_func == 'max':
                peaks = self.int_sorted_deisotoped_max_peaks
                ranks = self.deisotoped_max_ranks
            else:
                raise ValueError(f'Unknown deisotope_func: {deisotope_func}')
        else:
            peaks = self.intensity_sorted_peaks
            ranks = self.ranks

        if as_list:
            peak_idx = peaks.index(peak)
        else:
            peak_idx = peaks.index(peak.as_list())

        rank = ranks[peak_idx]
        return rank

    def calculate_sequence_coverage(self, lossy=False):
        """
        calculates the total, alpha, beta, bLike and yLike sequence coverage and the symmetry
        :return: dict
            'seq_coverage'
            'alpha_coverage'
            'beta_coverage'
            "symmetry"
            "alpha_pep_id"
            "beta_pep_id"
            "bLike_seq_coverage"
            "yLike_seq_coverage"
        """

        seq_cov_fragments = [f for f in self.get_fragments(lossy=lossy) if f.by_type != '']

        # alpha beta assignment
        pep1_fragments = [f for f in seq_cov_fragments if f.peptide_id == 0]
        pep2_fragments = [f for f in seq_cov_fragments if f.peptide_id == 1]

        n_pep1 = len(pep1_fragments)
        n_pep2 = len(pep2_fragments)
        from random import randrange
        # assign the alpha and beta peptide depending on which peptide has more fragments
        if n_pep1 > n_pep2:
            alpha = 0
            beta = 1
        elif n_pep2 > n_pep1:
            alpha = 1
            beta = 0
        # if the number of fragments is equal do random assignment of alpha and beta peptide
        elif randrange(2) == 0:
            alpha = 0
            beta = 1
        else:
            alpha = 1
            beta = 0

        aa_len_pep1 = sum(1 for c in self.peptide.pep_seq1 if c.isupper())
        aa_len_pep2 = sum(1 for c in self.peptide.pep_seq2 if c.isupper())

        seq_cov_fragments_set = set([f.sequence_coverage_id for f in seq_cov_fragments])
        seq_cov_frag_count = len(seq_cov_fragments_set)

        alpha_cov_fragments = set(
            [f.sequence_coverage_id for f in seq_cov_fragments if f.peptide_id == alpha])
        alpha_cov_frag_count = len(alpha_cov_fragments)

        beta_cov_fragments = set(
            [f.sequence_coverage_id for f in seq_cov_fragments if f.peptide_id == beta])
        beta_cov_frag_count = len(beta_cov_fragments)

        if alpha == 0:
            theoretical_alpha_frag_count = ((aa_len_pep1 - 1) * 2)
            theoretical_beta_frag_count = ((aa_len_pep2 - 1) * 2)
        else:
            theoretical_alpha_frag_count = ((aa_len_pep2 - 1) * 2)
            theoretical_beta_frag_count = ((aa_len_pep1 - 1) * 2)

        if self.isLinear:
            theoretical_frag_count = theoretical_alpha_frag_count
            beta_coverage = 0
        else:
            theoretical_frag_count = theoretical_alpha_frag_count + theoretical_beta_frag_count
            beta_coverage = float(beta_cov_frag_count) / theoretical_beta_frag_count

        seq_coverage = float(seq_cov_frag_count) / theoretical_frag_count
        alpha_coverage = float(alpha_cov_frag_count) / theoretical_alpha_frag_count

        if seq_cov_frag_count != alpha_cov_frag_count + beta_cov_frag_count:
            print("n: %s; n_a:%s; n_b: %s " % (seq_cov_frag_count, alpha_cov_frag_count,
                                               beta_cov_frag_count))
            raise Exception("Number of fragments doesn't add up")

        symmetry = 1 - math.fabs(alpha_coverage - beta_coverage)

        b_like_seq_cov_fragments = set(
            [f.sequence_coverage_id for f in seq_cov_fragments if f.by_type == 'bLike'])
        b_like_frag_count = len(b_like_seq_cov_fragments)

        y_like_seq_cov_fragments = set(
            [f.sequence_coverage_id for f in seq_cov_fragments if f.by_type == 'yLike'])
        y_like_frag_count = len(y_like_seq_cov_fragments)

        b_like_seq_cov = float(b_like_frag_count) / theoretical_frag_count
        y_like_seq_cov = float(y_like_frag_count) / theoretical_frag_count

        if b_like_frag_count + y_like_frag_count != seq_cov_frag_count:
            raise Exception("b/y like fragment counts don't add up!")

        # total_frag_count = len(found_fragments)
        # total_alpha_frag_count = len([f for f in found_fragments if f['pep'] == alpha])
        # total_beta_frag_count = len([f for f in found_fragments if f['pep'] == beta])
        # total_bLike_frag_count = len([f for f in found_fragments if f['type'] == 'bLike'])
        # total_yLike_frag_count = len([f for f in found_fragments if f['type'] == 'yLike'])

        return {
            'seq_coverage': seq_coverage,
            'alpha_coverage': alpha_coverage,
            'beta_coverage': beta_coverage,
            "symmetry": symmetry,
            # "n_frags_alpha": total_alpha_frag_count,
            # "n_frags_beta": total_beta_frag_count,
            # "n_frags": total_frag_count,
            "alpha_pep_id": alpha,
            "beta_pep_id": beta,
            # "dsso_alpha": dsso_alpha,
            # "dsso_beta": dsso_beta,
            "bLike_seq_coverage": b_like_seq_cov,
            "yLike_seq_coverage": y_like_seq_cov,
            # "n_frags_bLike": total_bLike_frag_count,
            # "n_frags_yLike": total_yLike_frag_count,
            # "alpha_LSdet": LSdet_alpha,
            # "beta_LSdet": LSdet_beta,
            # "ppm_errors": ppm_errors
        }


class Peptide:
    def __init__(
            self,
            pep_seq1,
            pep_seq2='',
            link_pos1=0,
            link_pos2=0,
    ):
        self.pep_seq1 = pep_seq1
        self.pep_seq2 = pep_seq2

        if pep_seq2 == '':
            self.is_linear = True
        else:
            self.is_linear = False

        self.link_pos1 = int(link_pos1)
        self.link_pos2 = int(link_pos2)

    def unique_id(self):
        if self.is_linear:
            return self.pep_seq1
        link_pos1 = str(self.link_pos1)
        link_pos2 = str(self.link_pos2)
        pep_seq1 = str(self.pep_seq1)
        pep_seq2 = str(self.pep_seq2)

        peptides = [
            {"seq": pep_seq1, "linkpos": link_pos1},
            {"seq": pep_seq2, "linkpos": link_pos2},
        ]
        peptides = sorted(peptides, key=lambda k: k['seq'])

        unique_id = "%s_%s-%s_%s" % (
            peptides[0]['seq'],
            peptides[1]['seq'],
            peptides[0]['linkpos'],
            peptides[1]['linkpos'],
        )

        return unique_id


class MzSpecies(Peptide):
    def __init__(
        self,
        pep_seq1,
        pep_seq2,
        charge,
        mz=None,
        link_pos1=0,
        link_pos2=0,
        rt=None
    ):
        Peptide.__init__(self, pep_seq1, pep_seq2, link_pos1, link_pos2)
        self.charge = charge
        # if mz is None:
        #     self.mz = self.calculate_mz()
        # else:
        #     self.mz = mz
        self.mz = mz
        self.rt = rt

    def unique_id(self):
        return "{}:{}".format(Peptide.unique_id(self), self.charge)

    # ToDo: write function
    # def calculate_mz(self):
    #     raise Exception('WIP')


class MzSpeciesTarget(MzSpecies):
    def __init__(
        self,
        pep_seq1,
        pep_seq2,
        charge,
        mz=None,
        link_pos1=0,
        link_pos2=0,
        rt_start=None,
        rt_stop=None
    ):
        MzSpecies.__init__(self, pep_seq1, pep_seq2, charge, mz, link_pos1, link_pos2)
        self.charge = int(charge)
        self.mz = float(mz)
        self.rt_start = float(rt_start)
        self.rt_stop = float(rt_stop)

    def match(self, mz, mz_tolerance, charge, rt, rt_delta=0):
        """
        matches the MzSpeciesTarget against a spectrum, ignores charge state requirement if
        charge is 0 (unassigned).
        :param mz: precursor mz
        :param mz_tolerance: mz matching tolerance
        :param charge: precursor charge state
        :param rt: spectrum rt
        :param rt_delta: rt error
        :return: True or False
        """
        if int(charge) != 0 and int(charge) != self.charge:
            return False

        rt_lo = self.rt_start - rt_delta
        rt_hi = self.rt_stop + rt_delta
        if not rt_lo <= rt <= rt_hi:
            return False

        mz_match_error = math.fabs((mz - self.mz) / self.mz)
        if mz_match_error > mz_tolerance:
            return False

        return True


class Peak:
    def __init__(self, p_id, mz, intensity, spectrum, cluster_ids=()):
        """
        Initialise the Peak.

        :param p_id: (int) peak id (index into the spectrum peak list)
        :param mz: (float) m/z value of the peak
        :param intensity: (float) intensity value of the peak
        :param spectrum: (Spectrum) spectrum object the peak belongs to
        :param cluster_ids: (list of int) cluster ids that the peak is part of
        """
        self.id = p_id
        self.mz = mz
        self.intensity = intensity
        self.spectrum = spectrum
        self.cluster_ids = cluster_ids

    def get_mz(self):
        return self.mz

    def get_intensity(self):
        return self.intensity

    def as_list(self):
        return [self.mz, self.intensity]

    def match_error(self, target_mz):
        return (target_mz - self.mz) / self.mz

    def get_rank(self):
        return self.spectrum.get_peak_rank(self)


class Fragment:
    def __init__(
            self,
            name,
            peptide_id,
            cluster,
            lossy,
            missing_monoisotopic,
            sequence,
            calc_mz,
            error,
            peak,
            spectrum,
            frag_type,
            charge
    ):
        self.name = name
        self.peptide_id = peptide_id
        self.lossy = lossy
        self.missing_monoisotopic = missing_monoisotopic
        self.cluster = cluster
        self.sequence = sequence
        self.calc_mz = float(calc_mz)
        self.error = error
        self.peak = peak
        self.spectrum = spectrum
        self.frag_type = frag_type
        self.charge = charge

    def get_intensity(self, deisotoped=False, deisotope_func='sum'):
        """
        Return the intensity of the peak or cluster.

        :param deisotoped: (bool) if True return the deisotoped intensity.
        :param deisotope_func: (str) function to use for deisotoped intensity.
            Valid values: - 'sum': summed up cluster intensity
                          - 'max': maximum peak intensity in the cluster
        :return:
        """
        if deisotoped:
            return self.cluster.get_intensity(deisotope_func)
        else:
            return self.peak.intensity

    def get_mz(self):
        return self.peak.mz

    def get_error_ppm(self):
        if self.error['unit'] == 'ppm':
            return self.error['value']
        else:
            return None

    def get_error_abs(self):
        if self.error['unit'] == 'Da':
            return self.error['value']
        else:
            return None

    @memoized_property
    def ion_type(self):
        # ToDo: might need rework based on fragment['type']
        if re.search('^[abcxyz]', self.name):
            return self.name[0]
        if re.search('P_dsso_[AST]', self.name):
            return 'dsso'
        if re.search('P_ucl_[ABCD]', self.name):
            return 'ucl'
        if self.name.startswith('P+P'):
            return 'Precursor'
        if self.name.startswith('P'):
            if self.spectrum.isLinear:
                return 'Precursor'
            else:
                return 'Peptide'

        raise Exception('Unknown fragment type in name: %s' % self.name)

    @memoized_property
    def ion_number(self):
        try:
            return re.match('[abcxyz]([0-9]+)(?:\\+P)?_?', self.name).groups()[0]
        except AttributeError:
            return ''

    @memoized_property
    def cluster_charge(self):
        return self.cluster.get_charge()

    @memoized_property
    def by_type(self):
        by_type = ''
        if re.search('^[abc]', self.name):
            by_type = 'bLike'
        elif re.search('^[xyz]', self.name):
            by_type = 'yLike'
        return by_type

    @memoized_property
    def sequence_coverage_id(self):
        """
        returns the unique identifier for calculating sequence coverage
        by_type-ion_number-peptide_id
        """
        return "-".join([self.by_type, self.ion_number, str(self.peptide_id)])

    def get_lossy(self):
        return self.lossy

    def get_rank(self):
        return self.peak.get_rank()

    def get_rel_int_base_peak(self, deisotoped=False):
        """
        get the relative intensity to the base peak of the spectrum
        :param deisotoped: deisotoping on/off
        :return: (float) relative intensity to base peak
        """

        peak_int = self.get_intensity(deisotoped=deisotoped)
        base_peak_int = self.spectrum.get_base_peak(deisotoped=deisotoped).intensity

        return peak_int / base_peak_int

    def get_rel_int_precursor(self, deisotoped=False, manual_match=False,
                              manual_match_tolerance=None):
        """
        Get the relative intensity to the precursor

        :param deisotoped: deisotoping on/off
        :param manual_match: False means only the matched precursor is considered,
            True tries to find the precursor manually
        :param manual_match_tolerance: tolerance for matching precursor.
            Default is spectrum fragment tolerance.
        :return: (float) relative intensity to precursor
        """

        peak_int = self.get_intensity(deisotoped=deisotoped)

        if manual_match:
            manual_precursor_match = self.spectrum.match_unfragmented_precursor_peak(
                deisotoped=deisotoped, tolerance=manual_match_tolerance)
            if manual_precursor_match[0] is None:
                precursor_int = 0
            else:
                precursor_int = manual_precursor_match[0].get_intensity()
        else:
            precursor_int = self.spectrum.get_unfragmented_precursor_intensity()

        try:
            return peak_int / precursor_int
        except ZeroDivisionError:
            return float('inf')

    def get_intensity_ratio(self, deisotoped=False):
        """
        Get the intensity ratio of intensity to all intensity.

        :param deisotoped: deisotoping on/off
        :return: float intensity ratio
        """
        # ToDo: total_spectrum_intensity deisotoped?
        # ToDo: sum([p.get_intensity() for p in self.spectrum.get_peaks(deisotoped=True)])
        # ToDo: due to duplicate clusters this is vastly higher than total_spectrum_intensity
        total_spectrum_intensity = sum([p.get_intensity() for p in self.spectrum.get_peaks()])
        return float(self.get_intensity(deisotoped=deisotoped)) / total_spectrum_intensity

    def as_dict(self):

        deisotoped_sum_rank = self.spectrum.get_peak_rank(
            [self.get_mz(), self.get_intensity(deisotoped=True, deisotope_func='sum')],
            deisotoped=True, deisotope_func='sum', as_list=True)
        deisotoped_max_rank = self.spectrum.get_peak_rank(
            [self.get_mz(), self.get_intensity(deisotoped=True, deisotope_func='max')],
            deisotoped=True, deisotope_func='max', as_list=True)
        return {
            'intensity': self.get_intensity(),
            'deisotoped_intensity_sum': self.get_intensity(deisotoped=True, deisotope_func='sum'),
            'deisotoped_intensity_max': self.get_intensity(deisotoped=True, deisotope_func='max'),
            'name': self.name,
            'pep_id': self.peptide_id,
            'calc_mz': self.calc_mz,
            'match_mz': self.get_mz(),
            'ppm': self.get_error_ppm(),
            'error_abs': self.get_error_abs(),
            'charge': self.charge,
            'cluster_charge': self.cluster_charge,
            'seq': self.sequence,
            'type': self.ion_type,
            'number': self.ion_number,
            'frag_type': self.frag_type,
            'by_type': self.by_type,
            'lossy': self.get_lossy(),
            'matched_missing_monoisotopic': self.missing_monoisotopic,
            'rank': self.get_rank(),
            'deisotoped_sum_rank': deisotoped_sum_rank,
            'deisotoped_max_rank': deisotoped_max_rank,
            'rel_int_base_peak': self.get_rel_int_base_peak(),
            # 'deisotoped_rel_int_base_peak': self.get_rel_int_base_peak(deisotoped=True),
            # 'rel_int_precursor': self.get_rel_int_precursor(),
            # 'deisotoped_rel_int_precursor': self.get_rel_int_precursor(deisotoped=True),
            # 'rel_int_precursor_manualMatch': self.get_rel_int_precursor(manual_match=True),
            # 'deisotoped_rel_int_precursor_manualMatch': self.get_rel_int_precursor(
            #     deisotoped=True, manual_match=True),
            # 'intensity_ratio': self.get_intensity_ratio(),
            # 'deisotoped_intensity_ratio': self.get_intensity_ratio(deisotoped=True),
        }


class IsotopeCluster:
    def __init__(self, c_id, peaks, charge):
        """
        Initialise the IsotopeCluster.

        :param c_id: (int) cluster id (index into the spectrum cluster list)
        :param peaks: (list of Peak) peaks associated with this cluster
        :param charge: (int) charge state assigned to the cluster
        """
        self.id = c_id
        self.peaks = peaks
        self.charge = charge

    @memoized_property
    def intensity_sum(self):
        return sum([p.intensity for p in self.peaks])

    @memoized_property
    def intensity_max(self):
        return max([p.intensity for p in self.peaks])

    def get_intensity(self, deisotope_func):
        if deisotope_func == 'sum':
            return self.intensity_sum
        elif deisotope_func == 'max':
            return self.intensity_max
        else:
            raise ValueError(f'Unknown deisotope_func: {deisotope_func}')

    def get_charge(self):
        return self.charge

    def get_peaks(self):
        return self.peaks

    def get_first_peak(self):
        return self.peaks[0]
