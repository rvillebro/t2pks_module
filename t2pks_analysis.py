# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the type II PKS module
    The analysis determines starter unit, chain length, molecular weight of the product and
    product class. It does so by identifying cyclases, ketoreductases, oxygenases,
    methylases and glycosyltransferases....
    If applicable the regioselectivity and function of identified auxilary and tailoring
    protiens will be predicted.

"""

import logging
import json
from typing import Any, Dict, List, Tuple, Set, Union  # used in comment type hints #pylint: disable=unused-import
from collections import defaultdict
from Bio.SearchIO import QueryResult, HSP

from antismash.config import ConfigType
from antismash.common import subprocessing
from antismash.common import path
from antismash.common import fasta
from antismash.common.hmmscan_refinement import refine_hmmscan_results, gather_by_query, HMMResult
from antismash.common.utils import get_hmm_lengths, get_fasta_lengths
from antismash.common.secmet import Record, Cluster, CDSFeature

from .results import T2PKS_Results, T2PKS_Prediction


def run_t2pks_hmmscan(cluster: Cluster) -> Dict[str, List[HMMResult]]:
    """ Runs hmmscan for type II PKS proteins on coding sequences in cluster

        Arguments:
            cluster: Cluster on which the type II PKS hmmscan shall be run

        Returns:
            a dictionary of key: cds and value: list of HMMResults, for hmmscan results of the cluster
    """
    cluster_fasta = fasta.get_fasta_from_features(cluster.cds_children)
    hmm_file = path.get_full_path(__file__, "data", "t2pks.hmm")
    hmm_results = subprocessing.run_hmmscan(hmm_file, cluster_fasta, opts=['--cut_tc'])
    hmm_lengths = get_hmm_lengths(hmm_file)
    return refine_hmmscan_results(hmm_results, hmm_lengths)

def run_starter_unit_blastp(cluster: Cluster, cds_hmm_hits: Dict[str, List[HMMResult]]) ->  Dict[str, List[HMMResult]]:
    """ Runs blastp on starter unit coding sequences in given cluster

        Arguments:
            cluster: Cluster on which the blastp shall be run
            cds_hmm_hits: HMMResults by cds from type II PKS hmmscan

        Returns:
            None if no starter unit cds are present otherwise a dictionary of key: cds and value: list of HMMresults, for blastp results of the cluster
    """
    starter_unit_cds = {}
    for cds, hmm_hits in cds_hmm_hits.items():
        starter_unit_hit_ids = [hit.hit_id for hit in hmm_hits if hit.hit_id in ['KSIII', 'AT', 'AMID', 'LIG']]
        if starter_unit_hit_ids:
            starter_unit_cds[cluster.parent_record.get_cds_by_name(cds)] = starter_unit_hit_ids

    if starter_unit_cds:
        blastp_results = []
        blastp_fasta_files = set()
        for cds, starter_unit_hit_ids in starter_unit_cds.items():
            query_sequence = fasta.get_fasta_from_features([cds])
            for hit_id in starter_unit_hit_ids:
                blast_database = path.get_full_path(__file__, 'data', hit_id)
                blastp_results.extend(subprocessing.run_blastp(blast_database, query_sequence))
                blastp_fasta_files.add(path.get_full_path(__file__, 'data', hit_id + '.fasta'))
        
        fasta_lengths = {}
        for fasta_file in blastp_fasta_files:
            fasta_lengths.update(get_fasta_lengths(fasta_file))

        return refine_hmmscan_results(blastp_results, fasta_lengths)
    
    return {}

def sum_predictions(predictions: List[Tuple[str, Union[str, None], float]], ret=None) -> Union[List[Tuple[str, int]], None]:
    """ Sum cds predictions by proposed function

        Arguments:
            cds_list: a list of T2PKS_CDSPredictions of which predictions wish to be summed

        Returns:
            None if cds_list is empty otherwose a sorted list tuples containing proposed function and summed bitscore by highest bitscore
    """
    if predictions:
        summed_bitscore = defaultdict(float)
        lowest_evalue = defaultdict(lambda: None)
        for ptype, pfunc, bitscore, evalue in predictions:
            if pfunc is not None:
                keys = pfunc.split('_')[-1].strip().split('/')
                for key in keys:
                    summed_bitscore[key] += bitscore/len(keys)
                    if lowest_evalue[key] is None:
                        lowest_evalue[key] = evalue
                    elif evalue < lowest_evalue[key]:
                        lowest_evalue[key] = evalue
        return sorted([(key, summed_bitscore[key], lowest_evalue[key]) for key in summed_bitscore], key=lambda tup: tup[1], reverse=True)

    if ret is not None:
        return ret
    
    return []

def predict_product_class(cds_predictions: Dict[str, List[Tuple[str, Union[str, None], float]]]) -> Set[str]:
    classification_file = path.get_full_path(__file__, "data", "classification.json")
    with open(classification_file, 'r') as file:
        classification = json.load(file)

    potential_classes = []
    for cds, predictions in cds_predictions.items():
        for ptype, pfunc, bitscore, evalue in predictions:
            if ptype in classification and pfunc != None:
                potential_classes.append(set(classification[ptype][pfunc]))

    if potential_classes:
        return set.intersection(*potential_classes)
    return set()

def predict_molecular_weight(results: T2PKS_Prediction) -> Union[List[int],None]:
    weights_file = path.get_full_path(__file__, "data", "weights.json")
    with open(weights_file, 'r') as file:
        weights = json.load(file)
    
    tailoring_mw = 0
    for cds, predictions in results.cds_predictions.items():
        for ptype, pfunc, bitscore, evalue in predictions:
            if ptype in list(weights.keys()):
                if ptype == 'CYC' and pfunc != None and '/' in pfunc:
                    tailoring_mw += 2 * weights[ptype]
                else:
                    tailoring_mw += weights[ptype] if ptype != 'HAL' else weights[ptype]['cl'] # weights[ptype][pfunc.split('_')[-1]]

    mws = {}
    if results.starter_unit:
        for unit, *_ in results.starter_unit:
            for elongations, *_ in results.malonyl_elongations:
                for n in elongations.split('|'):
                    mws[unit+'_'+n] = weights['starter_unit'][unit] + int(n) * weights['malonyl'] + tailoring_mw

    return mws

def analyze_t2pks_cluster(results: T2PKS_Prediction, cds_hmm_hits: Dict[str, List[HMMResult]], cds_blastp_hits: Dict[str, List[HMMResult]]):
    """ Analyze given type II PKS gene cluster based on HMM and blastp hits

        Arguments:
            results: T2PKS_Results instance to contain results
            cds_hmm_hits: type II PKS HMM hits by cds
            cds_blastp_hits: type II PKS blastp hits by cds
    """
    if cds_hmm_hits:
        for cds, hmm_hits in cds_hmm_hits.items():
            if cds in cds_blastp_hits.keys():
                hmm_hits.extend(cds_blastp_hits[cds])
            for hit in hmm_hits:
                info = hit.hit_id.split('_', 1)
                protein_type = info[0]
                protein_function = info[1] if len(info) == 2 else None
                results.cds_predictions[cds].append((protein_type, protein_function, hit.bitscore, hit.evalue))

    results.starter_unit = sum_predictions(results.get_predictions_by_protein_type(['KSIII', 'AT', 'AMID', 'LIG']), [('acetyl', 0, 0)])
    results.malonyl_elongations = sum_predictions(results.get_predictions_by_protein_type(['CLF']))
    results.product_class = predict_product_class(results.get_cds_predictions_by_predicted_protein_type(['CLF', 'CYC']))
    results.molecular_weight = predict_molecular_weight(results) if results.malonyl_elongations else {}

    

def t2pks_analysis(record: Record, results: T2PKS_Results, options: ConfigType) -> T2PKS_Results:
    """ Find all relvant genes for analysis

        Arguments:
            record: the record to search
            results: T2PKS_Results instance to contain results
            options: an antismash config object

        Returns:
            a single T2PKS_Results instance with analysis results for every type II PKS cluster
    """
    logging.info("Anlysing type II PKS clusters")
    t2pks_clusters = [cluster for cluster in record.get_clusters() if 't2pks' in cluster.products]
    if t2pks_clusters:
        for cluster in t2pks_clusters:
            t2pks_hmm_results = run_t2pks_hmmscan(cluster)
            t2pks_blastp_results = run_starter_unit_blastp(cluster, t2pks_hmm_results)
            results.cluster_predictions[cluster.get_cluster_number()] = T2PKS_Prediction()
            analyze_t2pks_cluster(results.cluster_predictions[cluster.get_cluster_number()], t2pks_hmm_results, t2pks_blastp_results)

        print(results)
        with open('res.json', 'w+') as json_file:
            json.dump(results.to_json(), json_file, indent='  ')
    else:
        logging.info("No type II PKS clusters to analyze")

    return results
