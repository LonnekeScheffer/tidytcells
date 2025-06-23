from tidytcells._standardized_gene_symbol.standardized_ig_symbol import (
    StandardizedIgSymbol,
)
from tidytcells._resources import VALID_MUSMUSCULUS_IG


class StandardizedMusMusculusIgSymbol(StandardizedIgSymbol):
    _synonym_dictionary = dict()
    _valid_ig_dictionary = VALID_MUSMUSCULUS_IG
