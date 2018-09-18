#!/bin/python
from typing import List
from .results import T2PKS_Results

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common import path
from antismash.common.layers import ClusterLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Cluster


def will_handle(products: List[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return bool(set(products).intersection({"t2pks"}))


def generate_sidepanel(cluster_layer: ClusterLayer, results: T2PKS_Results,
                       record_layer: RecordLayer, options_layer: OptionsLayer) -> str:
	""" Generate the sidepanel HTML with results from the type II PKS module """
	env = Environment(loader=FileSystemLoader(path.get_full_path(__file__, 'templates')),
					  autoescape=True, undefined=StrictUndefined)
	template = env.get_template('sidepanel.html')
	t2pks_layer = T2PKSLayer(results, cluster_layer.cluster_feature, record_layer)
	sidepanel = template.render(record=record_layer,
								cluster=t2pks_layer,
								results=results,
								options=options_layer)
	return sidepanel

class T2PKSLayer(ClusterLayer):

	def __init__(self, results: T2PKS_Results, cluster_feature: Cluster, record: RecordLayer) -> None:
		self.number = cluster_feature.get_cluster_number()
		self.record = record