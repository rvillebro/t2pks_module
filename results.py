# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

from collections import defaultdict
import logging
from typing import Any, Dict, List, Tuple, Set, Union  # pylint: disable=unused-import

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record, AntismashDomain, GeneFunction


class T2PKS_Prediction(object):
    __slots__ = ['starter_unit', 'malonyl_elongations', 'product_class', 'molecular_weight', 'cds_predictions']

    def __init__(self):
        self.starter_unit = [] # type: List[Tuple[str, float, float]]
        self.malonyl_elongations = [] # type: List[Tuple[str, float, float]]
        self.product_class = [] # type: Set[str]
        self.molecular_weight = {} # type: Dict[str, float]
        self.cds_predictions = defaultdict(list) # type: Dict[str, List[Tuple[str, Union[str, None], float, float]]]

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self) -> str:
        string = 'Starter unit: ' + str(self.starter_unit)
        string += '\nMalonyl elongations: ' + str(self.malonyl_elongations)
        string += '\nProduct class: ' + ','.join(self.product_class)
        string += '\nMolecular weight: ' + str(self.molecular_weight)
        
        string += '\nCDSs: {}'.format(len(self.cds_predictions))
        for cds, predictions in self.cds_predictions.items():
            string += '\n{}\n {}'.format(cds, '\n '.join(map(str, predictions)))
        
        return string

    def get_cds_predictions_by_predicted_protein_type(self, protein_types: List[str]) -> Dict[str, List[Tuple[str, Union[str, None], float]]]:
        ret = {}
        for cds, predictions in self.cds_predictions.items():
            cds_predicted_protein_types = [pred[0] for pred in predictions]
            if set(protein_types).intersection(set(cds_predicted_protein_types)):
                ret[cds] = predictions
        return ret

    def get_predictions_by_protein_type(self, protein_types: List[str]) -> Dict[str, List[Tuple[str, Union[str, None], float]]]:
        ret = []
        for cds, predictions in self.cds_predictions.items():
            for pred in predictions:
                if pred[0] in protein_types:
                    ret.append(pred)
        return ret

class T2PKS_Results(ModuleResults):
    """ The combined results of the type 2 PKS module """
    _schema_version = 1
    __slots__ = ["cluster_predictions"]

    def __init__(self, record_id: str):
        super().__init__(record_id)
        self.cluster_predictions = {} # type: Dict[int, T2PKS_Prediction]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        string = ''
        for cluster_id, prediction in self.cluster_predictions.items():
            string += 'Cluster {}\n'.format(cluster_id)
            string += str(prediction)
        return string

    def to_json(self) -> Dict[str, Any]:
        results = {"schema_version": self._schema_version,
                   "record_id": self.record_id,
                   "cluster_predictions": {}}

        for cluster_id, prediction in self.cluster_predictions.items():
            json_prediction = {}
            json_prediction['starter_unit'] = prediction.starter_unit
            json_prediction['malonyl_elongations'] = prediction.malonyl_elongations
            json_prediction['product_class'] = list(prediction.product_class)
            json_prediction['molecular_weight'] = prediction.molecular_weight
            json_prediction['cds_predictions'] = {}
            for cds, cds_predictions in prediction.cds_predictions.items():
                json_prediction['cds_predictions'][cds] = cds_predictions

            results['cluster_predictions'][cluster_id] = json_prediction

        return results

    @staticmethod
    def from_json(json: Dict[str, Any], _record: Record) -> "T2PKS_Results":
        assert "record_id" in json
        if json.get("schema_version") != T2PKS_Results._schema_version:
            logging.warning("Mismatching schema version, dropping results")
            return None
        results = T2PKS_Results(json["record_id"])

        json_cluster_predictions = json.get("cluster_predictions", {})
        for cluster_id, json_prediction in json_cluster_predictions:
            cluster_prediction = T2PKS_Prediction()
            for key, value in json_prediction.items():
                if key == 'starter_unit':
                    cluster_prediction.starter_unit = [tuple(pred) for pred in value]
                elif key == 'malonyl_elongations':
                    cluster_prediction.malonyl_elongations = [tuple(pred) for pred in value]
                elif key == 'product_class':
                    cluster_prediction.product_class = set(value)
                elif key == 'molecular_weight':
                    cluster_prediction.molecular_weight = value
                elif key == 'cds_predictions':
                    for cds, cds_predictions in value.items():
                        cluster_prediction.cds_predictions[cds] = [tuple(pred) for pred in cds_predictions]
            results[cluster_id] = cluster_prediction

        return results

    def add_to_record(self, record: Record) -> None:
        """ Save type II PKS prediction in record. Cluster prediction are saved in Cluster feature under qualifiers. Gene functions are added to CDS Feature
        """
        for cluster_number, t2pks_prediction in self.cluster_predictions.items():
            cluster = record.get_cluster(cluster_number)
            cluster._qualifiers['t2pks_starter_unit'] = ['{} (Score: {:.1f}; E-value: {:f})'.format(unit, bitscore, evalue) for unit, bitscore, evalue in t2pks_prediction.starter_unit]
            cluster._qualifiers['t2pks_malonyl_elongations'] = ['{} (Score: {:.1f}; E-value: {:f})'.format(clf, bitscore, evalue) for clf, bitscore, evalue in t2pks_prediction.malonyl_elongations]
            cluster._qualifiers['t2pks_product_class'] = list(t2pks_prediction.product_class)
            cluster._qualifiers['t2pks_molecular_wieght'] = ['{} (Da): {:.3f}'.format(comb, weight) for comb, weight in t2pks_prediction.molecular_weight.items()]

            for cds, predictions in t2pks_prediction.cds_predictions.items():
                for ptype, pfunc, bitscore, evalue in predictions:
                    cds_feature = record.get_cds_by_name(cds)
                    ann = ptype if pfunc is None else '{} {}'.format(ptype, pfunc)
                    cds_feature.gene_functions.add(GeneFunction.ADDITIONAL, 't2pks', '{} (Score: {:.1f}; E-value: {:f})'.format(ann, bitscore, evalue))