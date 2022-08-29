import pytest
from calphy.composition_transformation import CompositionTransformation

def test_composition():
	alc = CompositionTransformation("tests/conf1.data", {"Al": 500}, {"Al": 496, "Li": 1, "C": 2, "O": 1})
	assert alc.mappingdict['1-1'] == 1
	assert alc.mappingdict['1-2'] == 2
	assert alc.mappingdict['1-3'] == 3
	assert alc.mappingdict['1-4'] == 4
	assert len(alc.unique_mappings) == 4
	assert alc.unique_mapping_counts[0] == 496
	a, b = alc.update_pair_coeff("pair_coeff * * filename Al Li")
	assert a == 'pair_coeff * * filename Al Al Al Al'
	assert b == 'pair_coeff * * filename Al Li C O'