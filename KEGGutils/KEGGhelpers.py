# =============================================================================
# MISC HELPER FUNCTIONS
# =============================================================================

def push_backslash(stuff):
    """ push a backslash before a word, dumbest function ever"""
    
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

def shift_pos(pos, label_shift):
    """shift a pos by (sx, xy) pixels"""
    shiftx = label_shift[0]
    shifty = label_shift[1]
    pos2 = pos.copy()
    for key, position in pos2.items():
        pos2[key] = ( position[0] + shiftx, position[1] + shifty)
        
    return pos2

def shorten_labels(label_dict, n):
    """cmon does it really need a description"""
    shorten_dict = label_dict.copy()
    for key, item in label_dict.items():
        shorten_dict[key] = item[:n]
        
    return shorten_dict

