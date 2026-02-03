from typing import Optional
from tidytcells._utils.alignment import get_compatible_symbols
from tidytcells._resources import SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES


class MhGene:
    def __init__(self, original_input, error, gene_name=None, allele_designation=None, species=None):
        self._original_input = original_input
        self._error = error
        self._gene_name = gene_name
        self._allele_designation = allele_designation if allele_designation is not None and len(allele_designation) > 0 else None
        self._species = species

        self._highest_precision_symbol = self._gene_name

        if self._gene_name is not None and self._allele_designation is not None:
            self._highest_precision_symbol = f'{self._gene_name}*{":".join(self._allele_designation)}'

    def __str__(self):
        str_repr = self.symbol

        if str_repr is not None:
            return str_repr
        else:
            return ""

    @property
    def original_input(self) -> Optional[str]:
        return self._original_input

    @property
    def error(self) -> Optional[str]:
        return self._error

    @property
    def is_standardized(self) -> bool:
        return self.error is None

    @property
    def attempted_fix(self) -> Optional[str]:
        if not self.is_standardized:
            return self._highest_precision_symbol

    @property
    def symbol(self) -> Optional[str]:
        if self.is_standardized:
            return self._highest_precision_symbol

    @property
    def allele(self) -> Optional[str]:
        if self.is_standardized and self._allele_designation is not None and self._gene_name is not None:
            return f'{self._gene_name}*{":".join(self._allele_designation)}'

    @property
    def gene(self) -> Optional[str]:
        if self.is_standardized and self._gene_name is not None:
            return self._gene_name

    @property
    def species(self) -> str:
        return self._species


class HLAGene(MhGene):
    def __init__(self, original_input, error, gene_name=None, allele_designation=None):
        super().__init__(original_input, error, gene_name, allele_designation, species="homosapiens")

    @property
    def protein(self) -> Optional[str]:
        if self.is_standardized and self._allele_designation is not None and self._gene_name is not None:
            return f'{self._gene_name}*{":".join(self._allele_designation[:2])}'


class ReceptorGene:

    def __init__(self, original_input, error, gene_name=None, allele_designation=None, subgroup_name=None, species=None):
        self._original_input = original_input
        self._error = error
        self._gene_name = gene_name
        self._allele_designation = allele_designation
        self._subgroup_name = subgroup_name
        self._species = species

        self._highest_precision_symbol = None

        if self._gene_name is not None and self._allele_designation is not None:
            self._highest_precision_symbol = f"{self._gene_name}*{self._allele_designation}"
        elif self._gene_name is not None:
            self._highest_precision_symbol = self._gene_name
        elif self._subgroup_name is not None:
            self._highest_precision_symbol = self._subgroup_name

    def __str__(self):
        str_repr = self.symbol

        if str_repr is not None:
            return str_repr
        else:
            return ""

    @property
    def original_input(self) -> Optional[str]:
        return self._original_input

    @property
    def error(self) -> Optional[str]:
        return self._error

    @property
    def is_standardized(self) -> bool:
        return self.error is None

    @property
    def attempted_fix(self) -> Optional[str]:
        if not self.is_standardized:
            return self._highest_precision_symbol

    @property
    def symbol(self) -> Optional[str]:
        if self.is_standardized:
            return self._highest_precision_symbol

    @property
    def allele(self) -> Optional[str]:
        if self.is_standardized and self._allele_designation is not None and self._gene_name is not None:
            return f"{self._gene_name}*{self._allele_designation}"

    @property
    def gene(self) -> Optional[str]:
        if self.is_standardized:
            return self._gene_name

    @property
    def subgroup(self) -> Optional[str]:
        if self.is_standardized:
            return self._subgroup_name

    @property
    def locus(self) -> Optional[str]:
        if self.is_standardized:
            locus = self.symbol[0:3]
            if "/D" in self.symbol:
                locus += "/D"
            return locus

    @property
    def receptor_type(self):
        if self.is_standardized:
            return self.symbol[0:2]

    @property
    def gene_type(self) -> Optional[str]:
        if self.is_standardized:
            return self.symbol[3]

    @property
    def species(self) -> str:
        return self._species

    def get_all_alleles(self, enforce_functional=True):
        '''
        Get all alleles related to the standardized symbol

        :param enforce_functional: If 'true', only functional alleles are returned
        :return: A list of allele names
        '''
        if self.is_standardized:
            aa_dict = SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES[self.receptor_type][self.species]

            return get_compatible_symbols(self.symbol, aa_dict, self.gene_type, self.locus, enforce_functional)

    def get_aa_sequences(self, sequence_type="ALL", enforce_functional=True):
        '''
        Get amino acid sequence information related to the alleles of the standardized symbol

        :param sequence_type: which sequence to return. This can be:
                                For V genes: 'FR1', 'FR2', 'FR3', 'CDR1', 'CDR2', 'V-REGION'
                                For D genes: 'D-REGION'
                                For J genes: 'J-REGION', 'J-MOTIF'
                                Or 'ALL' to return all available sequences
        :param enforce_functional: If 'true', only information for functional alleles is returned
        :return: A dictionary with allele names as keys and sequences as values
                 When sequence_type is 'ALL', the result is a nested dictionary with allele names
                 as outer keys, sequence types as inner keys, and sequences as inner values.
        '''
        sequence_type = sequence_type.upper()
        sequence_type = sequence_type + "-IMGT" if sequence_type in {"FR1", "FR2", "FR3", "CDR1", "CDR2"} else sequence_type

        if self.is_standardized:
            aa_dict = SUPPORTED_RECEPTOR_SPECIES_AND_THEIR_AA_SEQUENCES[self.receptor_type][self.species]

            alleles_of_interest = self.get_all_alleles(enforce_functional)

            if sequence_type == "ALL":
                return {allele: aa_dict[allele] for allele in alleles_of_interest}
            else:
                return {allele: aa_dict[allele][sequence_type] if sequence_type in aa_dict[allele] else None
                        for allele in alleles_of_interest}


class Junction:

    def __init__(self, original_input, error, corrected_junction=None, species=None):
        self._original_input = original_input
        self._error = error
        self._corrected_junction = corrected_junction
        self._species = species

    def __str__(self):
        str_repr = self.junction

        if str_repr is not None:
            return str_repr
        else:
            return ""

    @property
    def original_input(self) -> Optional[str]:
        return self._original_input

    @property
    def error(self) -> Optional[str]:
        return self._error

    @property
    def is_standardized(self) -> bool:
        return self.error is None

    @property
    def attempted_fix(self) -> Optional[str]:
        if not self.is_standardized:
            return self._corrected_junction

    @property
    def junction(self) -> Optional[str]:
        if self.is_standardized:
            return self._corrected_junction

    @property
    def cdr3(self) -> Optional[str]:
        if self.is_standardized:
            if self._corrected_junction is not None and len(self._corrected_junction) > 2:
                return self._corrected_junction[1:-1]

    @property
    def species(self) -> str:
        return self._species
