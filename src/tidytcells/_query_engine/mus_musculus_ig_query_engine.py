from tidytcells._query_engine.tr_query_engine import TrIgQueryEngine
from tidytcells._resources import VALID_MUSMUSCULUS_IG


class MusMusculusIgQueryEngine(TrIgQueryEngine):
    _valid_gene_dictionary = VALID_MUSMUSCULUS_IG
