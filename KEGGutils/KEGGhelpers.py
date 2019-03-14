# =============================================================================
# MISC
# =============================================================================
def push_backslash(stuff):
    stuff_url = ""
    
    if stuff is None:
        stuff = ""
    else:
        stuff_url = "/" + stuff
        
    return stuff, stuff_url


def replace_dict_value(dictionary, old_value, new_value):
    """ Selectively replaces values in a dictionary
    
    Parameters:
        :dictionary(dict): input dictionary
        :old_value: value to be replaced
        :new_value: value to replace
        
    Returns:
        :output_dictionary (dict): dictionary with replaced values"""

    for key, value in dictionary.items():
        if value == old_value:
            dictionary[key] = new_value
    return dictionary