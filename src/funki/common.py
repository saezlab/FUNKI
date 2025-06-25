from functools import reduce

_colors = {
    'red': '#eb3f1b',
    'yellow': '#f6c510',
    'aqua': '#12cbd0',
    'teal': '#09a687',
    'blue': '#007fa0',
    'black': '#000000',
    'white': '#ffffff',
    'gray': '#8a8a8a',
    'lightgray': "#DADADA",
}


def is_numeric(var):
    '''
    Checks whether a value is numerical or not.
    '''

    if hasattr(var, '__iter__') and type(var) is not str:

        return all([is_numeric(i) for i in var])

    try:

        float(var)

        return True

    except ValueError:

        return False


def rget(dct, keys):
    '''
    Retrieves a value from nested dictionaries by diving recursively through a
    given list of keys. Returns the value if found, None otherwise
    '''

    def get(d, k):

        return d.get(k, {})

    return reduce(get, [dct] + keys) or None
