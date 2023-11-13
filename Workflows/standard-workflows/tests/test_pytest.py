import pytest
from standard_workflows import sc_funcs


# A dictionary with different parameters per dataset that shall inherit parameters from dict_1
dict_2 = {1: {'name': {'John': 'Jonny', 
                      'Mandoline': ('Mandy', 'ADD', 123)}, 
                'age': '27', 
                'sex': {
                        'Male': 'm', 
                        'Female': 'w', 
                        'Zwitter': 'z'}}, 
            2: {'name': 'Marie', 'age': '22', 'sex': 'Female'}}

# A dictionary with general parameters
dict_1 = {'project_params': {
    'name': {'Thomas': 'Tom', 
             'Mandoline': ('Mandy', 'Mandi', 555)},
    'sex': {'take_only': ['Male', 'Female'],
            'Female': 'woman'}
}}

def test_getpath():
    """Tests if we can get the paths to values or keys in a dictionary. Should work with simple types, lists and dicts.
    """
    assert sc_funcs.getpath(dict_2, 'Female', search_value = False) == (1, 'sex', 'Female')
    assert sc_funcs.getpath(dict_2, 'Female', search_value = True) == (2, 'sex')
    assert sc_funcs.getpath(dict_2, 'ADD', search_value = True) == (1, 'name', 'Mandoline')
    assert sc_funcs.getpath(dict_2, ('Mandy', 'ADD'), search_value = True) == (1, 'name', 'Mandoline')
    assert sc_funcs.getpath(dict_2, {'Male': 'm', 'Female': 'w', 'Zwitter': 'z'}, search_value = True) == (1, 'sex')
   
def test_mergepaths():
    """Tests if we can execute sth. like a right_join(dict1, dict2). 
    The 'ADD' keyword in dict_2 should be replaced with the values from dict_1 so that we get the sum of both dicts. Otherwise dict_2 should overwrite values from dict_1.
    The 'take_only' key in dict_2 should lead to subsetting the merged dict by the given values.
    """
    assert sc_funcs.merge_dicts(dict_1, dict_2[1]) == {'age': '27',
              'name': {'Thomas': 'Tom', 'John': 'Jonny', 'Mandoline': ['123', '555', 'Mandi', 'Mandy']},
              'sex': {'Female': 'w', 'Male': 'm'}, 'age': '27'}
