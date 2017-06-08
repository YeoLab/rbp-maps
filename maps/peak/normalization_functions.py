import math

def norm(some_list, num_events):
    """
    Normalizes each position by dividing by the number of events.
    
    Parameters
    ----------
    some_list
    num_events

    Returns
    -------

    """
    # some_list_ps = [x+1 for x in some_list] # remove pseudocount, uncomment to add back in
    normed_list = [float(x)/num_events for x in some_list]
    return normed_list

def std_error(some_list, num_events):
    """
    Returns standard deviation error given list of events.
    
    Parameters
    ----------
    some_list
    num_events

    Returns
    -------

    """
    devs = []
    p_list = norm(some_list, num_events)
    q_list = [1 - p for p in p_list]
    for p, q in zip(p_list, q_list):
        devs.append(dev(p, q, num_events))
    return devs

def dev(p, q, n):
    return math.sqrt(p*q) / math.sqrt(n)

def get_std_error_boundaries(hist, n):
    plus = [x + y for x, y in zip(norm(hist, n), std_error(hist, n))]
    minus = [x - y for x, y in zip(norm(hist, n), std_error(hist, n))]
    return plus, minus

